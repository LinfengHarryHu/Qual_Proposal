"""Fine-tune the selected MACE-MH-1 head without reconstructing its architecture."""

from collections import defaultdict
import json
from pathlib import Path
import sys

from ase.geometry import find_mic
from ase.io import read
from mace.data import AtomicData, Configuration
from mace.tools import AtomicNumberTable, torch_geometric
import numpy as np
import torch


HEAD = "oc20_usemppbe"
EPOCHS = 28
PATIENCE = 8
LEARNING_RATE = 2e-5


def projector(atoms):
    weights = np.zeros((len(atoms),3), dtype=np.float32)
    if atoms.info["reaction"] == "release":
        weights[[19, *range(120,156), *range(158,170)],2] = 1
    else:
        vector = find_mic(atoms.positions[84]-atoms.positions[132], atoms.cell, atoms.pbc)[0]
        weights[157] = vector/np.linalg.norm(vector)
    return torch.from_numpy(weights)


def make_graph(atoms, model, number_table):
    configuration = Configuration(
        atomic_numbers=atoms.numbers, positions=atoms.positions, cell=atoms.cell.array,
        pbc=tuple(atoms.pbc), properties={"forces":atoms.arrays["REF_forces"]},
        property_weights={"forces":1.0}, head=HEAD,
    )
    graph = AtomicData.from_config(configuration, number_table, model.r_max.item(), heads=model.heads)
    if graph.head.item() != model.heads.index(HEAD):
        raise ValueError("Graph selected an unintended pretrained head")
    graph.projection = projector(atoms)
    return graph


def loader(graphs, batch_size, shuffle=False):
    return torch_geometric.dataloader.DataLoader(graphs, batch_size=batch_size, shuffle=shuffle,
                                                 num_workers=0, drop_last=False)


def energy_offsets(model, atoms, graphs):
    model.eval()
    residuals = defaultdict(list)
    for structure, batch in zip(atoms, loader(graphs, 1)):
        with torch.no_grad():
            value = model(batch.to("cuda").to_dict(), compute_force=False)["energy"].item()
        residuals[structure.info["reaction"]].append(float(structure.info["REF_energy"])-value)
    return {group:float(np.mean(values)) for group,values in residuals.items()}


def errors(model, batch, training):
    batch = batch.to("cuda")
    pred = model(batch.to_dict(), training=training, compute_force=True)
    delta_energy = pred["energy"].reshape(-1)-batch.centered_energy.reshape(-1)
    delta_force = pred["forces"]-batch.forces
    directional_delta = torch.zeros(batch.ptr.numel()-1, device="cuda", dtype=delta_force.dtype)
    directional_delta.index_add_(0, batch.batch, (delta_force * batch.projection).sum(dim=1))
    return delta_energy, delta_force, directional_delta


def validate(model, samples, graphs):
    model.eval()
    energy, force, direction = [], [], []
    for batch in loader(graphs, 2):
        e, f, d = errors(model, batch, False)
        energy.extend(e.detach().tolist())
        force.append(f.detach().square().mean().item())
        direction.extend(d.detach().abs().tolist())
    result = {}
    for group in ("water","release"):
        pick = [i for i,a in enumerate(samples) if a.info["reaction"] == group]
        values = np.asarray(energy)[pick]
        result[group] = {
            "relative_energy_MAE_eV":float(np.mean(np.abs(values-values.mean()))),
            "directional_force_MAE_eV_A":float(np.mean(np.asarray(direction)[pick])),
        }
    score = sum(x["relative_energy_MAE_eV"]+x["directional_force_MAE_eV_A"]
                for x in result.values())/2.0
    return dict(groups=result, score=score, all_atom_force_RMSE_eV_A=float(np.sqrt(np.mean(force))))


def main(data, foundation, out):
    torch.manual_seed(42)
    torch.set_default_dtype(torch.float32)
    out.mkdir()
    model = torch.load(foundation, map_location="cpu", weights_only=False).float().to("cuda")
    if HEAD not in model.heads:
        raise ValueError("Checkpoint does not contain the intended head")
    numbers = AtomicNumberTable(model.atomic_numbers.tolist())
    samples = {split:read(data/f"{split}.extxyz", index=":") for split in ("train","valid")}
    assert {k:len(v) for k,v in samples.items()} == {"train":84,"valid":25}
    graphs = {k:[make_graph(atoms,model,numbers) for atoms in group] for k,group in samples.items()}
    offsets = energy_offsets(model,samples["train"],graphs["train"])
    for split in ("train","valid"):
        for atoms, graph in zip(samples[split], graphs[split]):
            value = float(atoms.info["REF_energy"])-offsets[atoms.info["reaction"]]
            graph.centered_energy = torch.tensor([value], dtype=torch.float32)
    print("Training-only energy offsets by composition (eV):", offsets, flush=True)
    baseline = validate(model,samples["valid"],graphs["valid"])
    print("validation baseline",json.dumps(baseline),flush=True)
    optimizer = torch.optim.AdamW(model.parameters(),lr=LEARNING_RATE,weight_decay=5e-7)
    history, best, best_epoch = [], float("inf"), 0
    for epoch in range(1,EPOCHS+1):
        model.train()
        losses = []
        for batch in loader(graphs["train"], 2, shuffle=True):
            optimizer.zero_grad(set_to_none=True)
            e, f, d = errors(model,batch,True)
            loss = 0.5*e.square().mean() + 3.0*f.square().mean() + 0.15*d.square().mean()
            if not torch.isfinite(loss):
                raise FloatingPointError(f"Nonfinite training loss in epoch {epoch}")
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(),5.0,error_if_nonfinite=True)
            optimizer.step()
            losses.append(loss.item())
        result = validate(model,samples["valid"],graphs["valid"])
        row = dict(epoch=epoch,training_loss=float(np.mean(losses)),**result)
        history.append(row)
        print(json.dumps(row,allow_nan=False),flush=True)
        if result["score"] < best:
            best, best_epoch = result["score"], epoch
            model.eval()
            torch.save(model,"best.model.tmp")
            Path("best.model.tmp").replace(out/"zno_mh1_pilot.model")
        if epoch-best_epoch >= PATIENCE:
            break
    if not best_epoch:
        raise ValueError("No validated fine-tuned checkpoint")
    meta = dict(head=HEAD, learning_rate=LEARNING_RATE, offsets_eV=offsets,
                validation_baseline=baseline, best_epoch=best_epoch,
                epochs_completed=len(history), history=history)
    (out/"validation.json").write_text(json.dumps(meta,indent=2,allow_nan=False)+"\n")


if __name__ == "__main__":
    data, foundation, target=sys.argv[1:]
    main(Path(data), foundation, Path(target))
