# PyInstaller recipe of the downloadable DNAKit desktop app; run it through build.py:
#
#   python desktop-app/build.py --install
#
# The result is dist/DNAKit/ (Windows, Linux) or dist/DNAKit.app (macOS). DNAKit must be
# installed as a regular (non-editable) package so that the bundled Skills are included.

import sys
from pathlib import Path

from PyInstaller.utils.hooks import collect_data_files, collect_submodules, copy_metadata

HERE = Path(SPECPATH)
ICON = str(HERE.parent / "docs" / "assets" / "images" / "DNAKit-icon.png")

# DNAKit loads most modules by name at run time and reads its own source and data files, so
# every submodule is included and the source files are kept next to the bytecode.
hiddenimports = collect_submodules("dnakit") + ["socksio"]
datas = collect_data_files("dnakit") + copy_metadata("dnakit")
for package in ("gradio", "gradio_client", "safehttpx", "groovy", "hf_gradio"):
    datas += collect_data_files(package)

a = Analysis(
    [str(HERE / "dnakit_desktop.py")],
    pathex=[],
    binaries=[],
    datas=datas,
    hiddenimports=hiddenimports,
    excludes=["tkinter", "pytest", "IPython", "matplotlib.tests", "numpy.tests", "pandas.tests"],
    module_collection_mode={"dnakit": "pyz+py", "gradio": "py", "gradio_client": "py"},
    noarchive=False,
)
pyz = PYZ(a.pure)
exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="DNAKit",
    console=sys.platform.startswith("linux"),
    icon=ICON,
    upx=False,
)
coll = COLLECT(exe, a.binaries, a.datas, name="DNAKit", upx=False)
if sys.platform == "darwin":
    app = BUNDLE(
        coll,
        name="DNAKit.app",
        icon=ICON,
        bundle_identifier="asia.mapengsen.dnakit",
        info_plist={"NSHighResolutionCapable": True, "CFBundleDisplayName": "DNAKit"},
    )
