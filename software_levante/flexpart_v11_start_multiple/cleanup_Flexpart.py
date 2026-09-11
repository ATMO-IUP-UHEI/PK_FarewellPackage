from pathlib import Path
import shutil
import datetime as dt
import pandas as pd


base_dir = Path("/work/bb1170/RUN/b383736/data/Flexpart_2021/Flexpart/insitu/")
start_date, end_date =dt.date(2019,10,1), dt.date(2023,3,31)

for month_start in pd.date_range(start_date, end_date, freq='MS'):
    config_dir=base_dir/str(month_start.strftime("%Y_%m"))
    print(config_dir)
    if not config_dir.is_dir():
        continue

    # Delete the entire config folder
    config_dir = config_dir / "config"
    if config_dir.exists():
        shutil.rmtree(config_dir)
        print(f"Deleted {config_dir}")

for day in pd.date_range(start_date,end_date):
    restart_dir=base_dir/str(day.strftime("%Y_%m"))/f'Release_{day.strftime("%Y%m%d")}'
    print(restart_dir)
    if restart_dir.exists():
        for file in restart_dir.glob("restart*"):
            if file.is_file():
                file.unlink()
                print(f"Deleted {file}")