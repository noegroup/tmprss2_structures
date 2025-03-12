import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.svm import SVR

morgan_df = pd.read_csv('morgan_fingerprints_drugbank.csv', header=None, index_col=0)
cddd_df = pd.read_csv('cddd_embeddings_drugbank.csv', header=None, index_col=0)
score_df = pd.read_csv('scores_md_drugbank.csv', index_col=0, usecols=[1, 2])
score_df.index = [idx[15:] for idx in score_df.index]
score_df = score_df.sort_index()
morgan_df = morgan_df[[idx in score_df.index for idx in morgan_df.index]]

# Parameters
db_size = score_df.shape[0]
n_clusters_percentage = 1 # percentage of dataset in initial set
ext_size_percentage = 1 # percentage of dataset in each extension
num_rounds = 9 # number of extension rounds
num_iterations = 5 # number of repetitions
score_type = "higher" # higher if higher score is better
random_selection = False # whether to select randomly or by training the model

n_clusters = round(db_size*n_clusters_percentage/100)
ext_size = round(db_size*ext_size_percentage/100)
np.random.seed(42)

for i in range(num_iterations):
    
    print('Iteration:', i+1, '\n')
    
    print("Clustering into "+str(n_clusters)+" clusters ("+str(n_clusters_percentage)+"%)")
    kmeans = KMeans(n_clusters=n_clusters, init="k-means++", n_init=20, max_iter=400, n_jobs=-1, algorithm="full")
    kmeans = kmeans.fit(morgan_df.values)
    cluster_df = pd.DataFrame(kmeans.labels_, morgan_df.index, columns=["cluster"])
    cluster_df.to_csv("clusters/clusters_"+str(n_clusters_percentage)+"_"+str(i)+".csv")

    # Get initial set
    current_ids = set(cluster_df['cluster'].drop_duplicates().index)
    potential_ids = set(morgan_df.index) - set(current_ids)

    for k in range(num_rounds):
        if random_selection:
            if len(potential_ids) >= ext_size:
                random_positions = np.random.choice(len(potential_ids), replace=False, size=ext_size)
            else:
                random_positions = np.arange(0, len(potential_ids))
            ext_ids = set(np.array(list(potential_ids))[random_positions])
        else:
            # get embeddings for current ids
            current_emb_df = cddd_df[[idx in current_ids for idx in cddd_df.index]]
            x_train = current_emb_df.values
            # get scores for current ids
            y_train = score_df[[idx in current_ids for idx in score_df.index]]['mscore'].values
            # get embeddings for potential ids
            potential_emb_df = cddd_df[[idx in potential_ids for idx in cddd_df.index]]
            x_pred = potential_emb_df.values
            # train and predict
            model = SVR(kernel='rbf', C=1.0)
            model.fit(x_train, y_train)
            y_pred = model.predict(x_pred)
            # update current and potential ids
            if score_type == "higher":
                ext_ids = set(potential_emb_df.iloc[np.argsort(y_pred)].index[-ext_size:]) # take last (higher score)
            else:
                ext_ids = set(potential_emb_df.iloc[np.argsort(y_pred)].index[:ext_size]) # take first
        
        current_ids = (current_ids | ext_ids)
        potential_ids = set(morgan_df.index) - set(current_ids)

    if score_type == "higher":
        score_sorted_df = score_df[[idx in current_ids for idx in score_df.index]].sort_values('mscore', ascending=False)
    else:
        score_sorted_df = score_df[[idx in current_ids for idx in score_df.index]].sort_values('mscore')
        
    score_sorted_df.to_csv("active_learning_round_"+str(i)+".csv")
