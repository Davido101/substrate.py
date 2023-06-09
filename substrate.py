import struct
import sys
from collections import OrderedDict
from io import BytesIO
import gzip

STREAM_MAGIC = bytes.fromhex('ac ed 00 05')
TC_BLOCKDATA = bytes.fromhex('77')

read_float = lambda f:struct.unpack('>f',f.read(4))[0]
read_int = lambda f:struct.unpack('>i',f.read(4))[0]
read_size = lambda f:ord(f.read(1))

def get_substrate_bytes(fname):
	data = b''

	with open(fname,'rb') as f:
		header = f.read(len(STREAM_MAGIC))

		if header != STREAM_MAGIC:
			print("Error: STREAM_MAGIC not found")
			exit()
		
		tc = f.read(1)

		if tc != TC_BLOCKDATA:
			print(f"Error: {tc} found instead of TC_BLOCKDATA")
			exit()
		else:
			data += f.read(read_size(f))

	return data

def read_struct(format,f):
    size = struct.calcsize(format)
    data = f.read(size)
    return list(struct.unpack(format,data))

def print_substrate_params(data):
	for key in data.keys():
		if key is None: continue
		print(f"{key}: {data[key]}")

def parse_substrate(data):
	substrate_fields = ('substrate_version','substrate_age','cell_count','environment_version','nutrient_rate','nutrient_chunk_size',
						'radiation_level','light_amount','light_direction_change','light_range',None,'cell_type_count',
						'spawn_phagocytes','spawn_flagellocytes','spawn_photocytes','spawn_devorocytes','spawn_lipocytes','spawn_keratinocytes',
						'spawn_buoyocytes','spawn_glueocytes','spawn_virocytes','spawn_nitrocytes','spawn_stereocytes','spawn_senseocytes',
						'spawn_myocytes','spawn_neurocytes','spawn_secrocytes','spawn_stemocytes','spawn_gametes','spawn_ciliocytes',
						'contaminate_with_random_cells','gravity','density','density_gradient','kill_cells_at_edge','nitrates',
						'max_cell_count','max_food_count','substrate_diameter','dynamic_friction','static_friction','only_point_mutations',
						'salinity','cell_aging','nutrient_lumpiness','nutrient_lump_size','mobile_food','nutrient_coating')

	data = BytesIO(data)
	version = read_int(data)
	if version != 95:
		print(f"Unsupported substrate version: {version}")

	params = [version]
	params += read_struct(">1d 2i 7d 1i",data)
	size = params[-1]
	params += read_struct(f">{size+1}? 3d 1? 1d 2i 3d 1? 1f 1? 2d 1? 1f",data)
	return OrderedDict(zip(substrate_fields,params))

substrate_data = get_substrate_bytes(sys.argv[1])
substrate_params = parse_substrate(substrate_data)
print_substrate_params(substrate_params)
with open(sys.argv[1],'rb') as f:
	f.seek(len(STREAM_MAGIC)+len(TC_BLOCKDATA)+1+len(substrate_data))
	compressed = f.read()

#print(gzip.decompress(compressed))