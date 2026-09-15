"""读 config.txt。所有脚本共用这一份。

config.txt 是 key = value 的纯文本，不引 YAML。路径一律按项目根目录解析，
所以从哪个目录调用脚本都能跑。

命令行也能用：python scripts/cfg.py workers    打印某个键的值，给 shell 脚本取参数用。
"""

from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CONFIG = ROOT / "config.txt"


def load(path=None):
    cfg = {}
    with open(path or CONFIG, encoding="utf-8") as fh:
        for raw in fh:
            line = raw.split("#", 1)[0].strip()
            if not line or "=" not in line:
                continue
            key, value = line.split("=", 1)
            cfg[key.strip()] = value.strip()
    return cfg


def i(cfg, key):
    return int(cfg[key])


def f(cfg, key):
    return float(cfg[key])


def lst(cfg, key):
    """逗号分隔的列表。"""
    return [x.strip() for x in cfg[key].split(",") if x.strip()]


def path(*parts):
    """项目内的相对路径。"""
    return ROOT.joinpath(*parts)


def out(*parts):
    """相对目录，不存在就建好。"""
    d = ROOT.joinpath(*parts)
    d.mkdir(parents=True, exist_ok=True)
    return d


def cjk_font():
    """让 matplotlib 能写中文。找不到中文字体就返回 None，图上改用英文标签。"""
    from matplotlib import font_manager, rcParams

    names = ["Noto Sans CJK SC", "Noto Sans CJK JP", "WenQuanYi Zen Hei",
             "Microsoft YaHei", "SimHei"]
    files = ["/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc",
             "/usr/share/fonts/truetype/wqy/wqy-zenhei.ttc",
             "/mnt/c/Windows/Fonts/msyh.ttc",
             "/mnt/c/Windows/Fonts/simhei.ttf"]
    have = {fo.name for fo in font_manager.fontManager.ttflist}
    for name in names:
        if name in have:
            rcParams["font.sans-serif"] = [name] + list(rcParams["font.sans-serif"])
            rcParams["axes.unicode_minus"] = False
            return name
    for fp in files:
        if Path(fp).exists():
            try:
                font_manager.fontManager.addfont(fp)
                name = font_manager.FontProperties(fname=fp).get_name()
                rcParams["font.sans-serif"] = [name] + list(rcParams["font.sans-serif"])
                rcParams["axes.unicode_minus"] = False
                return name
            except Exception:
                continue
    return None


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        for k, v in load().items():
            print(f"{k} = {v}")
    else:
        print(load()[sys.argv[1]])
