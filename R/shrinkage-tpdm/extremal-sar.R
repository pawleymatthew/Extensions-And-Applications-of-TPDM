# Function to create a binary rook neighborhood matrix
create_rook_matrix <- function(d) {
  # Total number of cells
  total_cells <- d * d
  
  # Create an empty matrix to store rook neighborhood relations
  rook_matrix <- matrix(0, nrow = total_cells, ncol = total_cells)
  
  # Define helper function to convert row, col index to a linear index
  index <- function(row, col) {
    return((row - 1) * d + col)
  }
  
  # Loop through each cell in the grid
  for (i in 1:d) {
    for (j in 1:d) {
      # Get the linear index for the current cell
      current_cell <- index(i, j)
      
      # Define rook neighbors: left, right, up, down
      if (i > 1) { # Up neighbor
        rook_matrix[current_cell, index(i - 1, j)] <- 1
      }
      if (i < d) { # Down neighbor
        rook_matrix[current_cell, index(i + 1, j)] <- 1
      }
      if (j > 1) { # Left neighbor
        rook_matrix[current_cell, index(i, j - 1)] <- 1
      }
      if (j < d) { # Right neighbor
        rook_matrix[current_cell, index(i, j + 1)] <- 1
      }
    }
  }
  
  return(rook_matrix)
}

make_dist_matrix <- function(grid_size) {
  coord <- expand.grid(seq(from = 0.5, to = grid_size - 0.5, by = 1),
                       seq(from = 0.5, to = grid_size - 0.5, by = 1))
  dist_matrix <- dist(coord, diag = TRUE, upper = TRUE) %>% as.matrix()
  return(dist_matrix)
}


sim_SAR <- function(n, A_tilde) {
  
  tau <- function(x, inverse = FALSE) {
    # softplus function and its inverse
    if (!inverse){
      y <- -plogis(x, lower.tail = FALSE, log.p = TRUE)
    } else {
      y <- qlogis(-x, lower.tail = FALSE, log.p = TRUE)
    }
    return(y)
  }
  
  d <- nrow(A_tilde)
  Z <- matrix(evd::rfrechet(n = n * d, shape = 2), nrow = d, ncol = n)
  X <- tau(A_tilde %*% tau(Z, inverse = TRUE)) %>% t()
  return(X)
}

lambda_soft_fix <- function(dist_matrix, Sigma_hat) {
  model <- nls(y ~ a * exp(-b * x) + c, 
               start = list(a = 1, b = 1, c = 0), # Initial parameter guesses
               data = data.frame(x = vecu(dist_matrix), y = vecu(Sigma_hat)))
  beta2 <- as.numeric(coefficients(model)["c"])
  return(beta2)
}

Sigma_rho <- function(rho, W) {
  d <- ncol(W)
  A_rho <- solve(diag(d) - rho * W)
  D_rho <- apply(A_rho, 1, function(x) sqrt(sum(x^2))) %>% diag()
  A_tilde_rho <- solve(D_rho) %*% A_rho
  Sigma_rho <- A_tilde_rho %*% t(A_tilde_rho)
  return(Sigma_rho)
}


