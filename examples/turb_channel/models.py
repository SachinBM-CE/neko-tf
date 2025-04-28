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
    def __init__(self, hidden_features=128):
        super(PolicyFunc, self).__init__()

        layers = [nn.Linear(in_features = 4,
                            out_features = hidden_features,
                            bias=True),
                  nn.ReLU(),
                  nn.Linear(in_features = hidden_features,
                            out_features = hidden_features // 2,
                            bias=True),
                  nn.ReLU(),
                  nn.Linear(in_features = hidden_features // 2,
                            out_features = 1,
                            bias=True),
                  nn.Tanh()]

        self.fwd = nn.Sequential(*layers)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.fwd(x)

class ValueFunc(nn.Module):
    def __init__(self, hidden_features=128):
        super(ValueFunc, self).__init__()

        layers = [nn.Linear(in_features = 5,
                            out_features = hidden_features,
                            bias=True),
                  nn.ReLU(),
                  nn.Linear(in_features = hidden_features,
                            out_features = hidden_features // 2,
                            bias=True),
                  nn.ReLU(),
                  nn.Linear(in_features = hidden_features // 2,
                            out_features = 1,
                            bias=True)]

        self.fwd = nn.Sequential(*layers)

    def forward(self, s: torch.Tensor, a: torch.Tensor) -> torch.Tensor:
        x = torch.cat([s, a], dim=1)
        return self.fwd(x)
