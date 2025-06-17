import math
import torch
from torch import nn
import torch.nn.functional as F

def weight_init(model, scale=0.02):
    with torch.no_grad():
        for m in model.modules():
            if isinstance(m, nn.Linear):
                sqrtk = math.sqrt(1./float(m.weight.shape[1]))
                nn.init.uniform_(m.weight, a=-sqrtk, b=sqrtk)
                if m.bias is not None:
                    m.bias.data.zero_()

class PolicyFunc(nn.Module):
    def __init__(self, hidden_features=128, action_dim=1): # action_dim defaults to 1 as confirmed
        super(PolicyFunc, self).__init__()

        # Define common layers for the policy network
        # Input features are 3 (state_dim) as per initialize_models.py
        self.common_layers = nn.Sequential(
            nn.Linear(in_features = 3, out_features = hidden_features, bias=True),
            nn.Softsign(),
            nn.Linear(in_features = hidden_features, out_features = hidden_features, bias=True),
            nn.Tanh()
        )

        # Output layer for the mean (mu) of the Gaussian distribution
        # The output dimension for mu is action_dim (which is 1)
        self.mu_layer = nn.Linear(in_features = hidden_features,
                                  out_features = action_dim,
                                  bias=True)

        # Output layer for the log standard deviation (log_sigma) of the Gaussian distribution
        # The output dimension for log_sigma is also action_dim (which is 1)
        self.log_sigma_layer = nn.Linear(in_features = hidden_features,
                                         out_features = action_dim,
                                         bias=True)

    def forward(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        # Pass the input through the common layers
        h = self.common_layers(x)
        
        # Get the mean (mu) from the mu_layer
        mu = self.mu_layer(h)
        
        # Get the log standard deviation (log_sigma) from the log_sigma_layer
        log_sigma = self.log_sigma_layer(h)
        log_sigma = torch.clamp(log_sigma, min=1e-6, max=1)

        # If you need to reshape to [batch_size, 1, 1]:
##        mu = mu.unsqueeze(-1)       # Adds a dimension at the end. Shape: [batch_size, 1, 1]
##        log_sigma = log_sigma.unsqueeze(-1) # Shape: [batch_size, 1, 1]

        # Return both mu and log_sigma as a tuple
        return mu, log_sigma

class ValueFunc(nn.Module):
    def __init__(self, hidden_features=128):
        super(ValueFunc, self).__init__()

        # Layers for the Value Function (Critic)
        # Input features are 4 (3 state_dim + 1 action_dim)
        layers = [nn.Linear(in_features = 4, out_features = hidden_features, bias=True),
                  nn.Softsign(),
                  nn.Linear(in_features = hidden_features, out_features = hidden_features, bias=True),
                  nn.Softsign(),
                  nn.Linear(in_features = hidden_features, out_features = 1, bias=True)]

        self.fwd = nn.Sequential(*layers)

    def forward(self, s: torch.Tensor, a: torch.Tensor) -> torch.Tensor:
        # Concatenate state and action tensors along dimension 1
        x = torch.cat([s, a], dim=1)
        # Pass the concatenated tensor through the sequential layers
        return self.fwd(x)

##import math
##import torch
##from torch import nn
##import torch.nn.functional as F
##
##def weight_init(model, scale=0.02):
##    with torch.no_grad():
##        for m in model.modules():
##            if isinstance(m, nn.Linear):
##                sqrtk = math.sqrt(1./float(m.weight.shape[1]))
##                nn.init.uniform_(m.weight, a=-sqrtk, b=sqrtk)
##                if m.bias is not None:
##                    m.bias.data.zero_()
##
##class PolicyFunc(nn.Module):
##    def __init__(self, hidden_features=128):
##        super(PolicyFunc, self).__init__()
##
##        layers = [nn.Linear(in_features = 3,
##                            out_features = hidden_features,
##                            bias=True),
##                  nn.Softsign(),
##                  nn.Linear(in_features = hidden_features,
##                            out_features = hidden_features,
##                            bias=True),
##                  nn.Softsign(),
##                  nn.Linear(in_features = hidden_features,
##                            out_features = 1,
##                            bias=True),
##                  nn.Softsign()]
##
##        self.fwd = nn.Sequential(*layers)
##
##    def forward(self, x: torch.Tensor) -> torch.Tensor:
####        x = torch.split(x, 3)
##        return self.fwd(x)
##
##class ValueFunc(nn.Module):
##    def __init__(self, hidden_features=128):
##        super(ValueFunc, self).__init__()
##
##        layers = [nn.Linear(in_features = 4,
##                            out_features = hidden_features,
##                            bias=True),
##                  nn.Softsign(),
##                  nn.Linear(in_features = hidden_features,
##                            out_features = hidden_features,
##                            bias=True),
##                  nn.Softsign(),
##                  nn.Linear(in_features = hidden_features,
##                            out_features = 1,
##                            bias=True),
##                  nn.Softsign()]
##
##        self.fwd = nn.Sequential(*layers)
##
##    def forward(self, s: torch.Tensor, a: torch.Tensor) -> torch.Tensor:
##        x = torch.cat([s, a], dim=1)
##        return self.fwd(x)
