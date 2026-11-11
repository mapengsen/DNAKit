Ready-to-run DNAKit desktop app. No Python is needed: download the file for your system, unpack it and start **DNAKit**.

| System | File | Start |
|---|---|---|
| Windows 10 / 11 (64-bit) | `DNAKit-windows-x64.zip` | Unpack, open the `DNAKit` folder, double-click `DNAKit.exe` |
| macOS, Apple silicon | `DNAKit-macos-arm64.dmg` | Open the disk image, drag **DNAKit** to Applications, open it |
| macOS, Intel | `DNAKit-macos-x64.dmg` | As above |
| Linux x86-64 (glibc 2.35 or later) | `DNAKit-linux-x64.tar.gz` | `tar -xzf DNAKit-linux-x64.tar.gz && ./DNAKit/DNAKit` |

The apps are not code-signed. At the first start, on Windows click *More info* → *Run anyway*; on macOS open *System Settings* → *Privacy & Security* and click *Open Anyway*.

Inputs and results are kept in the `DNAKit` folder of your home directory. The app contains the interface and DNAKit's core functions; functions that need optional packages (for example deep-learning models) need the pip installation. See the [graphical interface guide](https://github.com/mapengsen/DNAKit/blob/main/docs/gui.en.md).

These files are replaced whenever the desktop app is rebuilt from the `main` branch.
