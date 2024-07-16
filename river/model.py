import torch
import torch.nn as nn


class SimpleNN(nn.Module):
    def __init__(self, input_size, pos_dim, out_dim, hidden_dim=64):
        super(SimpleNN, self).__init__()
        self.exp_encoder = nn.Sequential(
            nn.Linear(input_size, hidden_dim),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(hidden_dim, hidden_dim)  # Binary classification
        )
        self.position_encoder =  nn.Sequential(
            nn.Linear(pos_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(hidden_dim, hidden_dim)  # Binary classification
        )
        self.classifier = nn.Sequential(
            nn.Linear(hidden_dim*2, hidden_dim),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(hidden_dim // 2, out_dim)  # Binary classification
        )
        
    def forward(self, x, position):
        x = self.exp_encoder(x)
        pos = self.position_encoder(position)
        x = torch.cat([x, pos], dim=-1)
        x = self.classifier(x)
        return x
