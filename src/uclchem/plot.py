import matplotlib.pyplot as plt


def plot_rate_summary(prod, dest, step, xlabel="rate", top_k_rates=5):
    fig, ax = plt.subplots(2, 1, sharex=True, figsize=(7, top_k_rates))
    prod.iloc[step].sort_values(ascending=False)[:top_k_rates].plot.barh(
        ax=ax[0], title="production", logx=True
    )
    dest.iloc[step].sort_values(ascending=False)[:top_k_rates].plot.barh(
        ax=ax[1], title="destruction"
    )
    ax[1].set_xlabel(xlabel)
