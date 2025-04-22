import yaml

def read_yaml_file(file_path):
    """
    Reads a YAML file and returns its content as a Python object (usually a dictionary or list).
    """
    try:
        with open(file_path, 'r') as file:
            return yaml.safe_load(file)
    except FileNotFoundError:
        print(f"Error: File not found: {file_path}. Please find example yaml in package directory.")
        return None
    except yaml.YAMLError as e:
         print(f"Error parsing YAML file: {e}")
         return None

file_path = 'config.yaml'
config = read_yaml_file(file_path)

# if config == None:
#     print("using default config")
#     config = {
#     "b": 3.9,
#     "c": 20,
#     "waveform": "smooth square",
#     "reflect": True,
#     "wave_amp": 0.7
#     }