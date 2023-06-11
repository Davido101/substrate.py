import struct
import sys
from collections import OrderedDict
from io import BytesIO
import gzip

STREAM_MAGIC = bytes.fromhex('ac ed 00 05')
TC_BLOCKDATA = bytes.fromhex('77')
TC_BLOCKDATALONG = bytes.fromhex('7a')

btd = lambda b:int(b.hex(),16)
read_float = lambda f:struct.unpack('>f',f.read(4))[0]
read_int = lambda f:struct.unpack('>i',f.read(4))[0]
read_double = lambda f:struct.unpack('>d',f.read(8))[0]
read_size = lambda f:btd(f.read(1))
read_size_long = lambda f:btd(f.read(4))

def error(msg):
	print(msg)
	exit()

def get_bytes(data,is_file=False,all_blocks=True):
	if not is_file: data = BytesIO(data)
	byt = b''

	header = data.read(len(STREAM_MAGIC))

	if header != STREAM_MAGIC:
		error("Error: STREAM_MAGIC not found")
	
	while True:
		tc = data.read(1)

		if not tc: break
		elif tc != TC_BLOCKDATA and tc != TC_BLOCKDATALONG:
			if not all_blocks: break
			error(f"Error: {tc} found instead of TC_BLOCKDATA")
		else:
			if tc == TC_BLOCKDATA:
				byt += data.read(read_size(data))
			else:
				byt += data.read(read_size_long(data))

	return byt

def get_file_bytes(fname,all_blocks=True):
	with open(fname,'rb') as f:
		return get_bytes(f,is_file=True,all_blocks=all_blocks)

def read_struct(form,f,big_endian=True):
	if big_endian:
		form = '>'+form

	size = struct.calcsize(form)
	data = f.read(size)
	return list(struct.unpack(form,data))

def zip_params(fields,params):
	return OrderedDict(zip(fields,params))

def print_params(data,depth=0):
	for key in data.keys():
		if key is None: continue
		factor = 1
		unit = ''
		value = data[key]
		if isinstance(key,tuple):
			factor = key[1]
			if len(key) >= 3:
				unit = key[2]

			key = key[0]

		key = key.replace('_',' ')
		if isinstance(value,bool): value = 'yes' if value else 'no'
		elif isinstance(value,list):
			print(f"{key}:")
			for e in value:
				print_params(e,depth+1)

			continue
		else: value *= factor

		prefix = '' if depth == 0 else '-'*depth+' '
		print(f"{prefix}{key}: {value}{unit}")

def parse_substrate(data):
	substrate_fields = ('substrate_version',('substrate_age',1,'h'),'cell_count','environment_version','nutrient_rate','nutrient_chunk_size',
						'radiation_level','light_amount','light_direction_change','light_range','substrate_viscosity','cell_type_count',
						'spawn_phagocytes','spawn_flagellocytes','spawn_photocytes','spawn_devorocytes','spawn_lipocytes','spawn_keratinocytes',
						'spawn_buoyocytes','spawn_glueocytes','spawn_virocytes','spawn_nitrocytes','spawn_stereocytes','spawn_senseocytes',
						'spawn_myocytes','spawn_neurocytes','spawn_secrocytes','spawn_stemocytes','spawn_gametes','spawn_ciliocytes',
						'contaminate_with_random_cells','gravity','density','density_gradient','kill_cells_at_edge','nitrates',
						'max_cell_count','max_food_count','substrate_diameter','dynamic_friction','static_friction','only_point_mutations',
						'salinity','cell_aging','nutrient_lumpiness','nutrient_lump_size','mobile_food','nutrient_coating')

	data = BytesIO(data)
	version = read_int(data)
	if version != 95:
		error(f"Unsupported substrate version: {version}")

	params = [version]
	params += read_struct("1d 2i 7d 1i",data)
	ntypes = params[-1]
	params += read_struct(f"{ntypes+1}? 3d 1? 1d 2i 3d 1? 1f 1? 2d 1? 1f",data)
	return zip_params(substrate_fields,params)

def parse_link(data):
	link_fields = (None,None,None,None,None,None,None,None)
	return zip_params(link_fields,read_struct("1i 2d 1? 2d 2f",data))

def parse_gene(data):
	gene_fields = (None,None,None,None,None,
				   None,None,None,None,None,
				   None,None,None,None,None,
				   None,None,None,None,None,
				   None,None,None,None,None,None)

	gene_params = read_struct("1i 9f 2i 3? 2i 4? 1f",data)
	extra0 = []
	for _ in range(12):
		extra0.append(read_struct("2h 3f",data))

	extra1 = read_struct("12i",data)
	gene_params += read_struct("2f",data)
	gene_params.append(extra0)
	gene_params.append(extra1)
	return zip_params(gene_fields,gene_params)

def parse_gzip(data,ncells):
	cell_fields = ('genome_version','x','y',None,None,
				   ('x_speed',500,'µm/h'),('y_speed',500,'µm/h'),None,None,('cell_diameter',1000,'µm'),
				   ('cell_mass',10,'ng'),('cell_age',1,'h'),'adhesin_connection_count','adhesin_connections',None,
				   None,None,None,None,'gene_count',
				   'genes',None,None,None,None,
				   None,None,('nitrogen_reserve',100,'%'),None,None,
				   None,None,None,None,None,
				   None,None,('toxins',100,'%'),('cell_injury',100,'%'),None,
				   None,None,('lipids',10,'ng'),None,None,None)

	data = BytesIO(data)
	param = read_double(data)
	cells = []

	for _ in range(ncells):
		genome_version = read_int(data)
		if genome_version != 95:
			error(f"Unsupported genome version: {genome_version}")

		cell_params = [genome_version]
		cell_params += read_struct("11d 1i",data)
		
		nlinks = cell_params[-1]
		links = []
		for __ in range(nlinks):
			links.append(parse_link(data))

		cell_params.append(links)

		cell_params += read_struct("1i 1? 3f 1i",data)
		ngenes = cell_params[-1]
		genes = []
		for __ in range(ngenes):
			genes.append(parse_gene(data))

		cell_params.append(genes)

		cell_params += read_struct("3i 4d 1? 14f 2i 1d",data)

		cells.append(zip_params(cell_fields,cell_params))


	return (param,cells)

if __name__ == '__main__':
	substrate_data = get_file_bytes(sys.argv[1],False)
	substrate_params = parse_substrate(substrate_data)
	print_params(substrate_params)
	with open(sys.argv[1],'rb') as f:
		f.seek(len(STREAM_MAGIC)+len(TC_BLOCKDATA)+1+len(substrate_data))
		compressed = f.read()

	with open('test.g','wb') as f:
		f.write(gzip.decompress(compressed))

	data = get_bytes(gzip.decompress(compressed))
	param,cells = parse_gzip(data,substrate_params['cell_count'])
	print()
	print(param)
	print()
	for cell in cells:
		input('Press Enter to see next cell')
		print()
		print_params(cell)
