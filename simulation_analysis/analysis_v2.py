import os


def list_directories(path):
    return [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]

def find_PT(string):
    if "PT" in string:
        return True
    else:
        return False
Size=18
devices = 1
for device in range(0, devices):
	path=f'N{Size}/DEV{device}'
	directories = list_directories(path)
	#print(directories)

	for d in directories:
		path_d = f'{path}/{d}'
		directory = list_directories(path_d)

		if find_PT(d):
			print(f'{path_d} is a Parallel Tempering Simulation')
		else:
			print(f'{path_d} is a Dynamical Simulation')
