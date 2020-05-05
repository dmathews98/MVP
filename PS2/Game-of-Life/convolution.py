k = np.array([[1,1,1],
			  [1,9,1],
	 	      [1,1,1]])
	rule = np.array([0,0,0,1,0,0,0,0,0,
					 0,0,1,1,0,0,0,0,0,0])
	return rule[signal.convolve2d(grid,k,boundary='wrap',mode='same').astype(int)]
