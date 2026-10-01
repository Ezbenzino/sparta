# ASCII only on purpose: Windows PowerShell 5.1 reads .ps1 as ANSI/GBK,
# so a non-BOM UTF-8 script with CJK text fails to parse.
# The real driver is Python. This file is just a thin wrapper.
& .\.venv\Scripts\python.exe .\scripts\run_fastpath.py $args
