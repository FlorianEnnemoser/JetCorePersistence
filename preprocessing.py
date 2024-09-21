from dataclasses import dataclass, field

@dataclass
class Region:
    name : str
    lat : tuple
    lon : tuple

    lat_tuple : tuple = field(init=False, default = None, repr=True)
    lon_tuple : tuple = field(init=False, default = None, repr=True)

    def __post_init__(self):
        self.lat_tuple : tuple = self.lat
        self.lon_tuple : tuple = self.lon
        self.lon_rect : list = [self.lon[0],self.lon[1],self.lon[1],self.lon[0],self.lon[0]]
        self.lat_rect : list = [self.lat[0],self.lat[0],self.lat[1],self.lat[1],self.lat[0]]
        self.lat : slice = slice(self.lat[0], self.lat[1])
        self.lon : slice = slice(self.lon[0], self.lon[1])