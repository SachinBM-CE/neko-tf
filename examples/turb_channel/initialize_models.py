import argparse as ap
import math
import torch
from functools import partial
from torch import nn
import torch.nn.functional as F

from models import weight_init, PolicyFunc, ValueFunc

def main(args):

    # set seed
    torch.manual_seed(666)

    # CUDA check
    if torch.cuda.is_available():
        torch.cuda.manual_seed(666)
        device = torch.device("cuda:0")
    else:
        device = torch.device("cpu")

    # parameters
    batch_size = 1164
    action_dim = 1 # Explicitly define action_dim as 1 here

    # policy model
    pmodel = PolicyFunc(hidden_features=args.num_hidden_features, action_dim=action_dim).to(device) # Pass action_dim
    weight_init(pmodel)
    jpmodel = torch.jit.script(pmodel)
    inp = torch.ones((batch_size, 3), dtype=torch.float32, device=device) # State input features: 3
    out_mu, out_log_sigma = jpmodel(inp) # Expect two outputs
    print("Policy model:", pmodel)
    print("Policy model mu output shape:", out_mu.shape)
    print("Policy model log_sigma output shape:", out_log_sigma.shape)
    torch.jit.save(jpmodel, "policy.pt")

    # value model
    qmodel = ValueFunc(hidden_features=args.num_hidden_features).to(device)
    weight_init(qmodel)
    jqmodel = torch.jit.script(qmodel)
    inp_a = torch.ones((batch_size, action_dim), dtype=torch.float32, device=device) # Action input features: 1
    out = jqmodel(inp, inp_a)
    print("Value model:", qmodel)
    print("Value model output shape:", out.shape)
    torch.jit.save(jqmodel, "value.pt")

if __name__ == "__main__":
    parser = ap.ArgumentParser()
    parser.add_argument("--num_hidden_features", type=int, default=128, help="Number of hidden features")
    args = parser.parse_args()

    main(args)

##import argparse as ap
##import math
##import torch
##from functools import partial
##from torch import nn
##import torch.nn.functional as F
##
##from models import weight_init, PolicyFunc, ValueFunc
##
##def main(args):
##
##    # set seed
##    torch.manual_seed(666)
##
##    # CUDA check
##    if torch.cuda.is_available():
##        torch.cuda.manual_seed(666)
##        device = torch.device("cuda:0")
##    else:
##        device = torch.device("cpu")
##
##    # parameters
##    batch_size = 1164
##
##    # policy model
##    pmodel = PolicyFunc(hidden_features=args.num_hidden_features).to(device)
##    weight_init(pmodel)
##    jpmodel = torch.jit.script(pmodel)
##    inp = torch.ones((batch_size, 3), dtype=torch.float32, device=device)
##    out = jpmodel(inp)
##    print("Policy model:", pmodel)
##    print("Policy model output shape:", out.shape)
##    torch.jit.save(jpmodel, "policy.pt")
##
##    # value model
##    qmodel = ValueFunc(hidden_features=args.num_hidden_features).to(device)
##    weight_init(qmodel)
##    jqmodel = torch.jit.script(qmodel)
##    inp_a = torch.ones((batch_size, 1), dtype=torch.float32, device=device)
##    out = jqmodel(inp, inp_a)
##    print("Value model:", qmodel)
##    print("Value model output shape:", out.shape)
##    torch.jit.save(jqmodel, "value.pt")
##
##if __name__ == "__main__":
##    parser = ap.ArgumentParser()
##    parser.add_argument("--num_hidden_features", type=int, default=128, help="Number of hidden features")
##    args = parser.parse_args()
##
##    main(args)
