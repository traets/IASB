nTrialFun <- function(nTripletsEff, distAttraction, distComp, distSimilarity, w1) {
  # author: nicolas skalicky, nicolas.skalicky@unibas.ch
  
  # this function creates an equal number of random attraction, compromise, and similarity sets
  # such that for half of the case A is the target (and B the competitor), for the other half - B is the target (and A the competitor)
  
  
  vecCoordinates = c() # the resulting vector of sample coordinates
  vecChoice = c() # the resulting vector of assumed choice preferences (1 = A, 2 = B, 3 = C)
  vecEffect = c() # the resulting vector of effects (1 = Attraction, 2 = Compromise, 3 = Similarity)
  
  w2 = 1 - w1 # attribute weight
  
  # PARAMETERS FOR THE RANDOM COORDINATES OF A 
  # (setting both to 0.5 lets A be randomised in the upper left quadrant)
  xmax = 0.5 # the maximal value of the x-coordinate at randomisation  --> x in the range (0 - 0.5)
  ymin = 0.5 # die minimal value of the y-coordinate at randomisation --> y in Range (0.5 - 1)
  
  # PARAMETERS FOR THE RANDOM COORDINATES OF B
  xmindiff = 0.1 # the minimal distance of B to A on the x-axis
  
  #go nTripletsEff-times through all different experiment settings and create coordinates
  
  # ATTRACTION
  if(nTripletsEff[1] != 0){
    for (i in 1:nTripletsEff[1]) {
      
      # put coordinates of C outside of area to enter while-loop
      Cx = 2
      Cy = 2
      
      while (Cx < 0 || Cx > 1 || Cy <0 || Cy > 1) {
        
        Ax <- runif(1, 0, xmax)
        Ay <- runif(1, ymin, 1)
        # calculate the y-distance at 0 of the straight line through A
        b = Ay + w1/w2*Ax
        
        # calculate the crossing of the straight line with the y-axis
        xzero = b*w2/w1
        
        # randomize B at any point in between A+the minimal distance and the crossing of the y-axis (or 1 if the crossing happens outside)
        Bx <- runif(1, Ax + xmindiff, min(xzero,1))
        By = -w1/w2*Bx + b
        
        # assure that every second attraction is set to point A
        # by assigning C to A for uneven numbers of i
        if (i %% 2 == 1) {
          # calculate x-coordinate o C and the y-axis distance at 0 for the orthogonal line
          Cx = Ax-(Bx-Ax)*distAttraction
          bortho = Ay-w2/w1*Ax
          
          # set preferred choice to A
          choice = 1
        }
        # assure that every second attraction is set to point B
        # by assigning C to B for even numbers of i
        if (i %% 2 == 0) {
          # calculate x-coordinate of C and the y-axis distance at 0 for the orthogonal line
          Cx = Bx-(Bx-Ax)*distAttraction
          bortho = By-w2/w1*Bx
          
          # set preferred choice to B
          choice = 2
        }
        
        
        Cy = w2/w1*Cx + bortho # calculate y-coordinate of C
        
      }
      
      # Add preferred choice to choice vector
      vecChoice = c(vecChoice, choice)
      # Add effect information to the corresponding vector: 1 = Attraction
      vecEffect = c(vecEffect, 1)
      # Add coordinates of A,B and C to the vector of all coordinates
      vecCoordinates = c(vecCoordinates, Ax, Bx, Cx, Ay, By, Cy)
      
    }
  }
  
  #COMPROMISE
  if(nTripletsEff[2] != 0){
    for (i in 1:nTripletsEff[2]) {
      
      # put coordinates of C outside of area to enter while-loop
      Cx = 2
      Cy = 2
      
      while (Cx < 0 || Cx > 1 || Cy <0 || Cy > 1) {
        Ax <- runif(1, 0, xmax)
        Ay <- runif(1, ymin, 1)
        # calculate the y-distance at 0 of the straight line through A
        b = Ay + w1/w2*Ax
        
        # calculate the crossing of the straight line with the y-axis
        xzero = b*w2/w1
        
        # randomize B at any point in between A+the minimal distance and the crossing of the y-axis (or 1 if the crossing happens outside)
        Bx <- runif(1, Ax + xmindiff, min(xzero,1))
        By = -w1/w2*Bx + b
        
        random <- runif(1, 0, 2)
        
        # assure that every second compromise is set to point A
        # by assigning C in the direction of A for uneven numbers of i
        if (i %% 2 == 1) {
          Cx = Ax-(Bx-Ax)*distComp
          # set preferred choice to A if A is in the middle (1=A)
          choice = 1
        }
        
        # assure that every second compromise is set to point B
        # by assigning C in the direction of B for even numbers of i
        if (i %% 2 == 0) {
          Cx = Bx+(Bx-Ax)*distComp
          # set preferred choice to B if C is in the middle (2=B)
          choice = 2
        }
        
        Cy = -w1/w2*Cx + b
        
      }
      
      # Add preferred choice to choice vector
      vecChoice = c(vecChoice, choice)
      # Add effect information to the corresponding vector: 3 = Similarity
      vecEffect = c(vecEffect, 2)
      # Add coordinates of A,B and C to the vector of all coordinates
      vecCoordinates = c(vecCoordinates, Ax, Bx, Cx, Ay, By, Cy)
      
    }
  }
  
  #SIMILARITY
  if(nTripletsEff[3] != 0){
    for (i in 1:nTripletsEff[3]) {
      
      
      # put coordinates of C outside of area to enter while-loop
      Cx = 2
      Cy = 2
      
      while (Cx < 0 || Cx > 1 || Cy <0 || Cy > 1) {
        Ax <- runif(1, 0, xmax)
        Ay <- runif(1, ymin, 1)
        # calculate the y-distance at 0 of the straight line through A
        b = Ay + w1/w2*Ax
        
        # calculate the crossing of the straight line with the y-axis
        xzero = b*w2/w1
        
        # randomize B at any point in between A+the minimal distance and the crossing of the y-axis (or 1 if the crossing happens outside)
        Bx <- runif(1, Ax + xmindiff, min(xzero,1))
        By = -w1/w2*Bx + b
        
        # assure that every second similiarity is set close to point A
        # by assigning C in the direction A for uneven numbers of i
        if (i %% 2 == 1) {
          Cx = Ax-(Bx-Ax)*distSimilarity
          # set preferred choice to B if C is close to A (2=B)
          choice = 2
        }
        # assure that every second similiarity is set close to point B
        # by assigning C in the direction B for even numbers of i
        if (i %% 2 == 0) {
          Cx = Bx+(Bx-Ax)*distSimilarity
          # set preferred choice to A if C is close to B (1=A)
          choice = 1
        }
        
        Cy = -w1/w2*Cx + b
        
      }
      
      # Add preferred choice to choice vector
      vecChoice = c(vecChoice, choice)
      # Add effect information to the corresponding vector: 2 = Compromise
      vecEffect = c(vecEffect, 3)
      # Add coordinates of A,B and C to the vector of all coordinates
      vecCoordinates = c(vecCoordinates, Ax, Bx, Cx, Ay, By, Cy)
      
    }
  }
  
  # merge all 3 vectors into the output list
  output = list(vecCoordinates, vecChoice, vecEffect)
  return(output)
  
} # end function: ntrial_fun