# Paramètres de simulation intrinsèques (sans angle injecté)
N = 60             # Nombre total de nœuds à générer
gamma = 0.05       # Taux d'expansion radial
m = 2              # Exposant du potentiel de répulsion
resolution = 1000  # Subdivision du cercle pour la recherche du minimum
import math

nodes = []  # Stockage des nœuds : { id, theta }

def computeEnergy(theta, t):
    energy = 0
    for k in range(len(nodes)):
        oldNode = nodes[k]
        r_k = math.exp(gamma * (t - oldNode['id']))
        cosDiff = math.cos(theta - oldNode['theta'])
        distSq = 1 + r_k * r_k - 2 * r_k * cosDiff
        if distSq > 1e-6:
            energy += 1 / (distSq ** (m / 2))
    return energy

def main():
    print("--- Début de la croissance auto-évitante ---")
    nodes.append({'id': 0, 'theta': 0})
    print("Nœud 00 -> Angle: 0.0000 | Divergence ρ: N/A")

    for t in range(1, N):
        minEnergy = float('inf')
        bestTheta = 0
        for i in range(resolution):
            candidateTheta = (i / resolution) * 2 * math.pi
            currentEnergy = computeEnergy(candidateTheta, t)
            if currentEnergy < minEnergy:
                minEnergy = currentEnergy
                bestTheta = candidateTheta

        prevTheta = nodes[-1]['theta']
        deltaTheta = (bestTheta - prevTheta) / (2 * math.pi)
        deltaTheta = deltaTheta - math.floor(deltaTheta)

        nodes.append({'id': t, 'theta': bestTheta})
        print(f"Nœud {str(t).zfill(2)} -> Angle: {bestTheta:.4f} | Divergence ρ: {deltaTheta:.6f}")

    finalRho = (nodes[N-1]['theta'] - nodes[N-2]['theta']) / (2 * math.pi)
    finalRho = finalRho - math.floor(finalRho)
    print("--------------------------------------------")
    print(f"Valeur asymptotique obtenue : ρ ≈ {finalRho:.6f}")
    print(f"Valeur théorique attendue (1/φ) : 0.618034 ou (1/φ²) : 0.381966")

if __name__ == '__main__':
    main()
