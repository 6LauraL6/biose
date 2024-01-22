import platform
import subprocess

def llista_directori():
    sistema_operatiu = platform.system()

    if sistema_operatiu == 'Windows':
        # Executa 'dir' en Windows
        subprocess.run(['dir'], shell=True)
    elif sistema_operatiu in ['Linux', 'Darwin']:  # 'Darwin' és el valor per a macOS
        # Executa 'ls -la' en Linux o macOS
        subprocess.run(['ls', '-la'])
    else:
        print(f"No es pot determinar el sistema operatiu: {sistema_operatiu}")

# Crida a la funció
llista_directori()
