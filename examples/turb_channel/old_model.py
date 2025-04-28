import torch
import torch.nn.functional as F

class Policy_Net(torch.nn.Module):
    def __init__(self, input_dim, action_dim, hidden_dim=128):
        """
        Neural network that outputs policy mean, policy standard deviation, and state value estimate.
        
        Args:
            input_dim (int) : Dimension of the input state.
            action_dim (int): Dimension of the action space (for the policy mean and std).
            hidden_dim (int): Number of units in each hidden layer.
        """
        super(Policy_Net, self).__init__()
        # Common backbone layers.
        self.fc1 = torch.nn.Linear(input_dim, hidden_dim)
        self.fc2 = torch.nn.Linear(hidden_dim, hidden_dim)
        
        # Policy heads.
        self.mean_head = torch.nn.Linear(hidden_dim, action_dim)
        self.std_head = torch.nn.Linear(hidden_dim, action_dim)
        
        # Value head.
        self.value_head = torch.nn.Linear(hidden_dim, 1)
    
    def forward(self, s):
        # First hidden layer with softsign activation.
        h1 = F.softsign(self.fc1(s))
        # Second hidden layer with softsign activation and a skip connection from h1.
        h2 = F.softsign(self.fc2(h1) + h1)
        
        # Policy outputs: mean and standard deviation.
        mu = self.mean_head(h2)
        # Use softplus to ensure standard deviation is positive.
        sigma = F.softplus(self.std_head(h2))
        
        # State value estimate.
        value = self.value_head(h2)
        
        return mu, sigma, value

def main():
    input_dim = 2    # Dimension of the input state.
    action_dim = 1   # Dimension of the action space.
    hidden_dim = 128 # Configurable hidden dimension.
    
    # Create the V-RACER network.
    model = Policy_Net(input_dim=input_dim, action_dim=action_dim, hidden_dim=hidden_dim)
    print("NN Model:", model)

    try:
        # Move model to GPU if available.
        model.to("cuda")
    except Exception as e:
        print("CUDA not available, running on CPU.")
    
    # Script and save the model.
    model_jit = torch.jit.script(model)
    model_jit.save("model.pt")

if __name__ == "__main__":
    main()
