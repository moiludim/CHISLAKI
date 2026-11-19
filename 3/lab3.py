import numpy as np
import matplotlib.pyplot as plt

class SmoothingSpline1D:
    
    def __init__(self, p=0.0):
        self.p = p              
        self.nodes_x = None    
        self.alpha = None       

    def update_spline(self, nodes_x, values_y, weights=None):
        self.nodes_x = np.asarray(nodes_x, dtype=float)
        y = np.asarray(values_y, dtype=float)
        n = len(self.nodes_x)
        
        if weights is None:
            weights = np.ones(n)
        else:
            weights = np.asarray(weights, dtype=float)

        A = np.zeros((n, n))
        b = np.zeros(n)

        term_data = 1.0 - self.p
        for i in range(n):
            A[i, i] += term_data * weights[i]
            b[i] += term_data * weights[i] * y[i]

        if self.p > 0.0:
            term_smooth = self.p
            for i in range(1, n):
                h_i = self.nodes_x[i] - self.nodes_x[i-1]
                if h_i <= 0:
                    continue
                inv_h = term_smooth / h_i
                
                A[i-1, i-1] += inv_h
                A[i-1, i]   -= inv_h
                A[i,   i-1] -= inv_h
                A[i,   i]   += inv_h

        self.alpha = np.linalg.solve(A, b)

    def get_value(self, x):
        return np.interp(x, self.nodes_x, self.alpha)


N_gen = 1963
M = 0.97
sigma = 5.24
nsmall = 200
PLOT_POINTS_PER_SEGMENT = 10

def main():
    try:
        np.random.seed(42)
        sample_y = np.random.normal(loc=M, scale=sigma, size=N_gen)
        sample_x = np.arange(1, N_gen + 1, dtype=float)

        if len(sample_x) < 2:
            raise ValueError("Недостаточно точек для сплайнов")

        count = min(nsmall, len(sample_x))
        nodes_x = sample_x[:count]
        sample_values = sample_y[:count]

        p_values = [0.0, 0.4, 0.8, 0.99]
        splines = {}

        for p in p_values:
            spline = SmoothingSpline1D(p=p)
            spline.update_spline(nodes_x, sample_values)
            splines[p] = spline

        plt.figure(figsize=(15, 8.5))

        plt.scatter(nodes_x, sample_values, color='black', s=20, label='Исходная выборка', zorder=5)

        x_dense = np.linspace(nodes_x[0], nodes_x[-1], (count - 1) * PLOT_POINTS_PER_SEGMENT + 1)

        colors = {0.0: 'black', 0.4: 'blue', 0.8: 'green', 0.99: 'orange'}
        
        for p in p_values:
            y_dense = np.array([splines[p].get_value(x) for x in x_dense])
            label_text = 'Интерполяционный сплайн (p=0)' if p == 0.0 else f'Сглаживающий сплайн p={p}'
            plt.plot(x_dense, y_dense, color=colors[p], linewidth=1.5, label=label_text)

        plt.title('Интерполяционный и сглаживающие сплайны (собственная реализация)', fontsize=11)
        plt.xlabel('Номер наблюдения', fontsize=9)
        plt.ylabel('Случайное значение', fontsize=9)
        plt.xlim(1, 50)
        plt.ylim(-40, 40)
        plt.grid(True)
        plt.legend(loc='upper right', bbox_to_anchor=(1.25, 1), fontsize=8)
        plt.tight_layout()
        plt.show()

    except Exception as e:
        print(f"Ошибка: {e}")

if __name__ == "__main__":
    main()