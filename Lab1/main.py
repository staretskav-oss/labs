import requests
import numpy as np
import matplotlib.pyplot as plt

# =====================================================
# 1. Запит до API висот
# =====================================================

url = "https://api.open-elevation.com/api/v1/lookup?locations=" \
      "48.164214,24.536044|48.164983,24.534836|48.165605,24.534068|" \
      "48.166228,24.532915|48.166777,24.531927|48.167326,24.530884|" \
      "48.167011,24.530061|48.166053,24.528039|48.166655,24.526064|" \
      "48.166497,24.523574|48.166128,24.520214|48.165416,24.517170|" \
      "48.164546,24.514640|48.163412,24.512980|48.162331,24.511715|" \
      "48.162015,24.509462|48.162147,24.506932|48.161751,24.504244|" \
      "48.161197,24.501793|48.160580,24.500537|48.160250,24.500106"

response = requests.get(url)
data = response.json()

results = data["results"]
print("Кількість вузлів:", len(results))

# =====================================================
# 2. Табуляція
# =====================================================

print("\n№ | Latitude | Longitude | Elevation")

for i, p in enumerate(results):
    print(i, "|", p["latitude"], "|", p["longitude"], "|", p["elevation"])

# =====================================================
# 3. Запис у файл
# =====================================================

with open("nodes.txt", "w") as f:
    f.write("№ | Latitude | Longitude | Elevation\n")
    for i, p in enumerate(results):
        f.write(f"{i} | {p['latitude']} | {p['longitude']} | {p['elevation']}\n")

# =====================================================
# 4. Кумулятивна відстань
# =====================================================

def haversine(lat1, lon1, lat2, lon2):
    R = 6371000
    phi1 = np.radians(lat1)
    phi2 = np.radians(lat2)

    dphi = np.radians(lat2 - lat1)
    dlambda = np.radians(lon2 - lon1)

    a = np.sin(dphi/2)**2 + np.cos(phi1)*np.cos(phi2)*np.sin(dlambda/2)**2
    return 2 * R * np.arctan2(np.sqrt(a), np.sqrt(1-a))


coords = [(p["latitude"], p["longitude"]) for p in results]
elevations = [p["elevation"] for p in results]

distances = [0]

for i in range(1, len(coords)):
    d = haversine(coords[i-1][0], coords[i-1][1], coords[i][0], coords[i][1])
    distances.append(distances[-1] + d)

distances = np.array(distances)
elevations = np.array(elevations)

print("\nDistance | Elevation")
for i in range(len(distances)):
    print(f"{distances[i]:.2f} | {elevations[i]:.2f}")

# =====================================================
# 5. Графік маршруту
# =====================================================

plt.figure(figsize=(10,6))
plt.plot(distances, elevations, marker='o', label="GPS точки")

plt.xlabel("Distance (m)")
plt.ylabel("Elevation (m)")
plt.title("Маршрут Заросляк -> Говерла")

plt.grid(True)
plt.legend()

plt.show()

# =====================================================
# 6. Кубічний сплайн
# =====================================================

def cubic_spline(x,y):

    n = len(x)
    h = np.diff(x)

    A = np.zeros((n,n))
    B = np.zeros(n)

    A[0,0] = 1
    A[n-1,n-1] = 1

    for i in range(1,n-1):
        A[i,i-1] = h[i-1]
        A[i,i] = 2*(h[i-1] + h[i])
        A[i,i+1] = h[i]

        B[i] = 3*((y[i+1]-y[i])/h[i] - (y[i]-y[i-1])/h[i-1])

    c = np.linalg.solve(A,B)

    a = y[:-1]
    b = np.zeros(n-1)
    d = np.zeros(n-1)

    for i in range(n-1):
        b[i] = (y[i+1]-y[i])/h[i] - h[i]*(2*c[i] + c[i+1])/3
        d[i] = (c[i+1]-c[i])/(3*h[i])

    return a,b,c[:-1],d


a,b,c,d = cubic_spline(distances,elevations)

# =====================================================
# 7. Побудова гладкого сплайна
# =====================================================

xx = np.linspace(distances[0], distances[-1], 500)
yy = []

for val in xx:

    for i in range(len(distances)-1):

        if distances[i] <= val <= distances[i+1]:

            dx = val - distances[i]
            y = a[i] + b[i]*dx + c[i]*dx**2 + d[i]*dx**3

            yy.append(y)
            break

yy = np.array(yy)

plt.figure(figsize=(10,6))

plt.plot(distances, elevations, 'o', label="GPS точки")
plt.plot(xx, yy, label="Cubic spline")

plt.xlabel("Distance (m)")
plt.ylabel("Elevation (m)")
plt.title("Spline Route Profile")

plt.grid(True)
plt.legend()

plt.show()

# =====================================================
# 8. Сплайн 10 / 15 / 20 вузлів
# =====================================================

plt.figure(figsize=(10,6))

for m in [10,15,20]:

    x_nodes = np.linspace(distances[0], distances[-1], m)
    y_nodes = np.interp(x_nodes, distances, elevations)

    a2,b2,c2,d2 = cubic_spline(x_nodes,y_nodes)

    xx2 = np.linspace(x_nodes[0], x_nodes[-1], 400)
    yy2 = []

    for val in xx2:

        for i in range(len(x_nodes)-1):

            if x_nodes[i] <= val <= x_nodes[i+1]:

                dx = val - x_nodes[i]
                y = a2[i] + b2[i]*dx + c2[i]*dx**2 + d2[i]*dx**3
                yy2.append(y)
                break

    plt.plot(xx2, yy2, label=f"{m} вузлів")

plt.scatter(distances, elevations, color="black", label="Original data")

plt.xlabel("Distance (m)")
plt.ylabel("Elevation (m)")
plt.title("Spline comparison (10 / 15 / 20 nodes)")

plt.grid(True)
plt.legend()

plt.show()

# =====================================================
# ДОДАТКОВЕ ЗАВДАННЯ
# =====================================================

print("\n=== Характеристики маршруту ===")

print("Загальна довжина маршруту (м):", distances[-1])

total_ascent = sum(max(elevations[i]-elevations[i-1],0) for i in range(1,len(elevations)))
print("Сумарний набір висоти (м):", total_ascent)

total_descent = sum(max(elevations[i-1]-elevations[i],0) for i in range(1,len(elevations)))
print("Сумарний спуск (м):", total_descent)

print("Мінімальна висота (м):", min(elevations))
print("Максимальна висота (м):", max(elevations))

# =====================================================
# Аналіз градієнта
# =====================================================

print("\n=== Аналіз градієнта ===")

grad = np.gradient(yy, xx) * 100

print("Максимальний підйом (%):", np.max(grad))
print("Максимальний спуск (%):", np.min(grad))
print("Середній градієнт (%):", np.mean(np.abs(grad)))

steep = xx[np.abs(grad) > 15]

print("Кількість точок з крутизною >15%:", len(steep))

if len(steep) > 0:
    print("Перша точка крутизни >15%:", round(steep[0],2), "м")
    print("Остання точка крутизни >15%:", round(steep[-1],2), "м")

# =====================================================
# Енергія підйому
# =====================================================

print("\n=== Механічна енергія підйому ===")

mass = 80
g = 9.81

energy = mass * g * total_ascent

print("Механічна робота (Дж):", energy)
print("Механічна робота (кДж):", energy/1000)
print("Енергія (ккал):", energy/4184)

