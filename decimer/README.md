## install PySide6
```powershell
 & 'C:\Program Files\Cresset-BMD\Flare\pyflare.exe' -m pip install --user pyside6
```

## How to Use Decimer in Flare Python Extensions

1. Run the following command on the server side:
```
python decimer_api.py
```

2. Place the file decimer_client_flare_ext.py into the Python extensions directory of Flare.

##  extracts chemical structure depictions from a PDF document
```
# 基本用法（页码从1开始，化合物编号从1开始）
python decimer-segmentation.py -i input.pdf -o compounds_screenshots

# 指定起始页码和起始化合物编号
python decimer-segmentation.py -i input.pdf -o output --start-page 40 --start-compound 101 --dpi 300
```
