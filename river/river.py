
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader, TensorDataset
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.preprocessing import LabelEncoder
from tqdm import tqdm
from captum.attr import IntegratedGradients, DeepLiftShap, GradientShap
from tqdm import  tqdm
from .model import SimpleNN




class River:
    def __init__(self, gene_expression, spatial, label):
        self.gene_expression = gene_expression
        self.spatial = spatial
        self.label = label
        self.device = 'cuda'
        self.model = SimpleNN(gene_expression.shape[1], spatial.shape[1], np.unique(label).shape[0]).to(self.device)
        self.final_rank = None 
        self.dataset = TensorDataset(torch.from_numpy(gene_expression).float(), torch.from_numpy(spatial).float(), torch.from_numpy(label))
        self.train_loader = DataLoader(self.dataset, batch_size=4096, shuffle=True)
        self.test_loader = DataLoader(self.dataset, batch_size=4096, shuffle=True)
        self.optimizer = optim.Adam(self.model.parameters(), lr=1e-3, weight_decay=1e-4)
        self.loss = nn.CrossEntropyLoss()
        self.n_features = gene_expression.shape[1]
        self.scores_ig = None
        self.scores_dl = None
        self.scores_sl = None
    
    def train(self, epoch=200):   
        for epoch in range(epoch):
            for data, pos, targets in self.train_loader:
                outputs = self.model(data.to(self.device), pos.to(self.device))
                loss = self.loss(outputs, targets.to(self.device))
                # Backward and optimize
                self.optimizer.zero_grad()
                loss.backward()
                self.optimizer.step()

            print(f'Epoch [{epoch+1}/100], Loss: {loss.item():.4f}')

            # Test the model
            self.model.eval()
            with torch.no_grad():
                correct = 0
                total = 0
                for data, pos, targets in self.test_loader:
                    outputs = self.model(data.to(self.device), pos.to(self.device))
                    _, predicted = torch.max(outputs.data, 1)
                    #print(predicted)
                    total += targets.size(0)
                    correct += (predicted.cpu() == targets.cpu()).sum().item()

            print('Accuracy of the network on the test data: %d %%' % (100 * correct / total))
        
    def attribution(self):
        torch.set_num_threads(32)
        ig = IntegratedGradients(self.model)
        dl = DeepLiftShap(self.model)
        sl = GradientShap(self.model)

        ig_attributions = []
        dl_attributions = []
        sl_attributions = []

        test_loader = DataLoader(self.dataset, batch_size=32, shuffle=False, drop_last=True)

        self.model.eval()

        self.model.to('cuda:3')

        for data, position, label in tqdm(test_loader):
            
            ig_attribution = ig.attribute((data.to('cuda:3'), position.to('cuda:3')), baselines=(torch.zeros(32, self.n_features).to('cuda:3'), position.to('cuda:3')), target=label.to('cuda:3'))
            
            #ig_attributions += torch.nn.functional.normalize(abs(ig_attribution.cpu()).reshape(-1, adata.X.shape[1]), dim=-1).sum(0)
            ig_attributions.append(abs(ig_attribution[0].cpu().detach()))
            
            del ig_attribution

            dl_attribution = dl.attribute((data.to('cuda:3'), position.to('cuda:3')), baselines=(torch.zeros(32, self.n_features).to('cuda:3'), position.to('cuda:3')), target=label.to('cuda:3'))
            
            #dl_attributions += torch.nn.functional.normalize(abs(dl_attribution.cpu()).reshape(-1, adata.X.shape[1]), dim=-1).sum(0)
            dl_attributions.append(abs(dl_attribution[0].cpu().detach()))

            del dl_attribution

            sl_attribution = sl.attribute((data.to('cuda:3'), position.to('cuda:3')), baselines=(torch.zeros(32, self.n_features).to('cuda:3'), position.to('cuda:3')), target=label.to('cuda:3'))
            
            sl_attributions.append(abs(sl_attribution[0].cpu().detach()))
            
            del sl_attribution
            
            torch.cuda.empty_cache()    
        
        return ig_attributions, dl_attributions, sl_attributions
        
    def summary_attribution(self, ig_attributions, dl_attributions, sl_attributions, overlap):
        final_score_ig =  torch.cat(ig_attributions).reshape(-1, self.n_features)
        final_score_dl =  torch.cat(dl_attributions).reshape(-1, self.n_features)
        final_score_sl =  torch.cat(sl_attributions).reshape(-1, self.n_features)
        scores_ig = torch.nn.functional.normalize(final_score_ig, dim=1).detach()
        scores_dl = torch.nn.functional.normalize(final_score_dl, dim=1).detach()
        scores_sl = torch.nn.functional.normalize(final_score_sl, dim=1).detach()

        ig_rank_indices = torch.argsort(-scores_ig.sum(0))
        dl_rank_indices = torch.argsort(-scores_dl.sum(0))
        sl_rank_indices = torch.argsort(-scores_sl.sum(0))

        print(np.array(overlap)[ig_rank_indices][:10].tolist())
        print(np.array(overlap)[dl_rank_indices][:10].tolist())
        print(np.array(overlap)[sl_rank_indices][:10].tolist())

    
        ig_rank = torch.argsort(ig_rank_indices)
        dl_rank = torch.argsort(dl_rank_indices)
        sl_rank = torch.argsort(sl_rank_indices)   
        
        data = {
            'Rank1': ig_rank.tolist(),
            'Rank2': sl_rank.tolist(),
            'Rank3': dl_rank.tolist()
        }
        df = pd.DataFrame(data)

        df.index = overlap#adata.var_names.values##adata.var_names#adata.var['label']

        # The lower the rank, the higher the score (e.g., 1st place gets highest score)
        n = len(df)
        borda_scores = df.apply(lambda x: n - x).sum(axis=1)

        # Sorting individuals by Borda scores in descending order
        final_ranks = borda_scores.sort_values(ascending=False)

        self.final_rank = final_ranks
        self.scores_ig = scores_ig
        self.scores_dl = scores_dl
        self.scores_sl = scores_sl
        
        
    def return_top_k_gene(self, top_k=200):
        return self.final_rank[:200].index.values 
        
    def save_model(self, path):
        torch.save(path, self.model)
        
    def load_model(self, path):
        self.model.load_state_dict(torch.load(path)['state_dict'])
        
    

