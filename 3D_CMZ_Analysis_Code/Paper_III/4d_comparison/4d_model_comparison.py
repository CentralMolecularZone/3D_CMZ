import numpy as np
import pandas as pd
from collections import Counter
from scipy.spatial.distance import cdist
from sklearn.preprocessing import RobustScaler
from sklearn.neighbors import NearestNeighbors


def load_data(path, sep='\t', names=['l', 'b', 'v', 'near_far']):
    return pd.read_csv(path, sep=sep, header=None, names=names)


def preprocess_data(df):
    df = df.copy()
    df['near_far_numeric'] = df['near_far'].map({'Near': 0, 'Far': 1})
    scaler = RobustScaler()
    normalized = pd.DataFrame(
        scaler.fit_transform(df[['l', 'b', 'v']]),
        columns=['l', 'b', 'v']
    )
    normalized['near_far'] = df['near_far']
    normalized['near_far_numeric'] = df['near_far_numeric']
    return normalized, df


def load_and_preprocess_models():
    model_files = [
        ('molinari_resampled_300.txt', '\t', "Molinari"),
        ('sofue_resampled_300.txt', '\t', "Sofue"),
        ('kdl_resampled_300.txt', '\t', "KDL"),
        ('ellipse_resampled_300.txt', '\t', "Ellipse")
    ]
    return [(*preprocess_data(load_data(file, sep)), name) for file, sep, name in model_files]


def calculate_mahalanobis_distances(data1, data2):
    cov = np.cov(data1.T)
    inv_cov = np.linalg.inv(cov)
    distances = cdist(data1, data2, metric='mahalanobis', VI=inv_cov)
    return distances, inv_cov


def predict_near_far(model_data, catalogue_data, n_neighbors, inv_cov):
    nn = NearestNeighbors(n_neighbors=n_neighbors, metric='mahalanobis', metric_params={'VI': inv_cov})
    nn.fit(model_data[['l', 'b', 'v']].values)
    distances, indices = nn.kneighbors(catalogue_data)

    predicted_nf = []
    weights = []
    for idx, dist in zip(indices, distances):
        neighbors_nf = model_data['near_far'].iloc[idx].values
        predicted_nf.append(Counter(neighbors_nf).most_common(1)[0][0])
        weight = 1 / np.median(dist)
        weights.append(weight)

    return np.array(predicted_nf), np.array(weights)


def analyse_model_4d_with_k(model, catalogue, model_name, k_range):
    model_data = model[['l', 'b', 'v']].values
    catalogue_data = catalogue[['l', 'b', 'v']].values

    distances, inv_cov = calculate_mahalanobis_distances(model_data, catalogue_data)
    normalized_distances = distances / np.max(distances)
    overall_distance_3d = np.mean(np.min(normalized_distances, axis=0))

    nf_accuracies = []
    for k in k_range:
        predicted_nf, weights = predict_near_far(model, catalogue_data, k, inv_cov)
        actual_nf = catalogue['near_far'].values
        nf_accuracy = np.sum((predicted_nf == actual_nf) * weights) / np.sum(weights)
        nf_accuracies.append(nf_accuracy)

    return {
        'name': model_name,
        'overall_distance_3d': overall_distance_3d,
        'nf_accuracies': nf_accuracies,
        'k_range': list(k_range),
    }


def run_knn(model, catalogue, model_name, k):
    model_data = model[['l', 'b', 'v']].values
    catalogue_data = catalogue[['l', 'b', 'v']].values

    distances, inv_cov = calculate_mahalanobis_distances(model_data, catalogue_data)
    predicted_nf, weights = predict_near_far(model, catalogue_data, k, inv_cov)
    actual_nf = catalogue['near_far'].values
    nf_accuracy = (np.sum((predicted_nf == actual_nf) * weights) / np.sum(weights)) * 100

    return {
        'name': model_name,
        'k': k,
        'nf_accuracy': nf_accuracy
    }


def main():
    catalogue, catalogue_original = preprocess_data(load_data('walker-catalogue.txt', sep=','))
    models = load_and_preprocess_models()

    k = int(np.sqrt(len(models[0][0])))
    knn_results = [run_knn(model, catalogue, name, k) for model, original_data, name in models]

    print("\nModel ranking:")
    sorted_results = sorted(knn_results, key=lambda x: x['nf_accuracy'], reverse=True)
    for i, result in enumerate(sorted_results, 1):
        print(f"{i}. {result['name']}: Near/Far accuracy = {result['nf_accuracy']:.0f}%")


if __name__ == "__main__":
    main()