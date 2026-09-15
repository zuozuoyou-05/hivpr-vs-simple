"""单个配体的对接。batchdock.sh 会把它并行调起来。

用法
  python scripts/dock_one.py --make-maps                      算一次格点
  python scripts/dock_one.py --ligand <in.pdbqt> --out <out.pdbqt> --log <out.log>

格点只跟受体和盒子有关，跟配体无关，所以整批配体共用一套（--make-maps 跑一次就够），
每个进程直接 load_maps，省掉每个配体重算格点那部分开销。

打分写成一份简单的两列表，extract_scores.py 按列解析。表头尽量贴近 vina 自己的
log 格式，方便肉眼对照；只是没有抄 rmsd 那两列，因为 Python API 的 energies()
里带出来的 rmsd 是负值（内部比较用的中间量，不是到最优构象的距离），写了反而是误导。
"""

import argparse
import json
import time
from pathlib import Path

import cfg

LOG_HEAD = """\
#################################################################
# vina 打分（hivpr-vs-simple）
#################################################################

 mode |   affinity
      |  (kcal/mol)
-----+------------
"""


def make_maps(receptor, box_path, prefix, cpu):
    from vina import Vina

    box = json.loads(Path(box_path).read_text(encoding="utf-8"))
    prefix = Path(prefix)
    # Vina 不肯覆盖已有的格点文件，先清干净
    for old in prefix.parent.glob(prefix.name + ".*"):
        if old.suffix == ".map" or old.name.endswith(".maps"):
            old.unlink()
    v = Vina(sf_name="vina", cpu=cpu, seed=0, verbosity=0)
    v.set_receptor(str(receptor))
    v.compute_vina_maps(center=box["center"], box_size=box["size"], force_even_voxels=True)
    v.write_maps(str(prefix))
    return box


def dock(args):
    from vina import Vina

    t0 = time.time()
    v = Vina(sf_name="vina", cpu=args.cpu, seed=args.seed, verbosity=0)
    v.load_maps(args.maps)
    v.set_ligand_from_file(args.ligand)
    v.dock(exhaustiveness=args.exh, n_poses=args.poses)
    energies = v.energies(n_poses=args.poses)

    Path(args.out).write_text(v.poses(n_poses=1), encoding="utf-8")
    rows = [LOG_HEAD]
    for k, row in enumerate(energies, 1):
        rows.append(f"{k:4d}    {row[0]:9.3f}  {row[1]:9.3f}  {row[2]:9.3f}\n")
    rows.append(f"# seconds {time.time() - t0:.1f}\n")
    Path(args.log).write_text("".join(rows), encoding="utf-8")
    return float(energies[0][0])


def main():
    cfgd = cfg.load()
    ap = argparse.ArgumentParser()
    ap.add_argument("--make-maps", action="store_true", help="只算格点不跑配体")
    ap.add_argument("--ligand")
    ap.add_argument("--out")
    ap.add_argument("--log")
    ap.add_argument("--maps", default=str(cfg.out("vina_outs") / "receptor"))
    ap.add_argument("--receptor", default=str(cfg.path("data", "processed", "receptor.pdbqt")))
    ap.add_argument("--box", default=str(cfg.path("data", "processed", "box.json")))
    ap.add_argument("--exh", type=int, default=cfg.i(cfgd, "exhaustiveness"))
    ap.add_argument("--poses", type=int, default=cfg.i(cfgd, "n_poses"))
    ap.add_argument("--seed", type=int, default=cfg.i(cfgd, "vina_seed"))
    ap.add_argument("--cpu", type=int, default=1)
    args = ap.parse_args()

    if args.make_maps:
        box = make_maps(args.receptor, args.box, args.maps, cpu=4)
        print(f"格点写出 {(args.maps)}.*，center={box['center']} size={box['size']}")
    else:
        aff = dock(args)
        print(f"{Path(args.ligand).stem}  {aff:.2f} kcal/mol")


if __name__ == "__main__":
    main()
