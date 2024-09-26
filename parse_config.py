import tomllib
from preprocessing import Region

def load_config(filepath : str) -> dict:
    with open(filepath, "rb") as f:
        config = tomllib.load(f)
    return config

def execute(filepath : str):
    config = load_config(filepath)

    region = Region(name=config["REGION"]["name"],
                    lat=config["REGION"]["latitude"],
                    lon=config["REGION"]["longitude"])
    
