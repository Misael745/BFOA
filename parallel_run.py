import time
import csv
from parallel_BFOA import run_bfoa

def experimento_bfoa(n_corridas=30, archivo_salida="resultados.csv"):
    with open(archivo_salida, "w", newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Corrida", "Fitness", "Tiempo (s)", "Interacciones", "BLOSUM"])

        for i in range(n_corridas):
            print(f"\n✨ Ejecutando corrida {i+1} de {n_corridas}...")
            start = time.time()
            fitness, interacciones, blosum = run_bfoa(iteraciones=30)
            tiempo = time.time() - start
            writer.writerow([i+1, fitness, tiempo, interacciones, blosum])
            print(f"✅ Corrida {i+1}: Fitness={fitness}, Tiempo={tiempo:.2f}s, Interacciones={interacciones}, BLOSUM={blosum}")

if __name__ == "__main__":
    experimento_bfoa()
