import os

def criar_arquivo_config(configpath):
    file = "config.yaml"
    create = True
    if os.path.exists(os.path.join(configpath, file)):
        print(f"A config file already exists at {os.path.join(configpath, file)}, you wish to overwrite it? (y/[n])")
        ow = input().strip().lower()
        if ow != 'y':
            print("Exiting without changes.")
            create = False
            return
        else:
            print("Overwriting existing config file.")
    if create:
        with open(os.path.join(configpath, file), 'w') as f:
            f.write("organisms:\n")
            f.write("k:\n")
            f.write("imp: false")
        print(f"Config file created in: {configpath}")

def main():
    os.makedirs('workflow/config', exist_ok=True)
    configpath = 'workflow/config'
    criar_arquivo_config(configpath)


if __name__ == "__main__":
    main()