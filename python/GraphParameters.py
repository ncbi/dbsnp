class GraphParameter:
    """
    Set parameters for the scatterplot
    
    """

    def __init__(self, args):
        score_error = self.ValidScores(args)
        if score_error:
            print(f'\n{score_error}\n')
            exit()
            
        self.xcol = args.x_score  # score to be plotted on x-axis
        self.ycol = args.y_score  # score to be plotted on y-axis

        self.width = 8            # width of graph in inches
        self.height = 8           # height of graph in inches

        self.lgd_x_pos = 0.02     # legend x position (0 = left, 0.5 = middle)
        self.lgd_y_pos = 0.03     # legend y position (0 = top, 0.5 = middle)
        self.lgd_size = 1.0       # legend size 
        self.lgd_font_size = 8    # legend font size
        self.lgd_x_gap = 1.5      # space between legend and label
        self.lgd_y_gap = 1.5      # space between two legend rows
        self.tick_lbl_size = 8    # tick label font size
        
        # Vertices of the EAF (European, African, East Asian) triangle
        self.e_x, self.e_y = 1.4748, 1.4370
        self.f_x, self.f_y = 1.0800, 1.1000
        self.a_x, self.a_y = 1.7045, 1.1000
        
        min_size, max_size = 3, 100
        if args.gw:
            gw = int(args.gw)
            if gw >= min_size and gw <= max_size:
                self.width = gw
                
        if args.gh:
            gh = int(args.gh)
            if gh >= min_size and gh <= max_size:
                self.height = gh

        if args.lgd_x:
            self.lgd_x_pos = float(args.lgd_x)
        if args.lgd_y:
            self.lgd_y_pos = float(args.lgd_y)

        ranges = self.GetDefaultAxisRanges()
        
        xrange = ranges[self.xcol]
        yrange = ranges[self.ycol]
        
        self.xmin = float(xrange[0])  # x-axis min val to be plotted
        self.xmax = float(xrange[1])  # x-axis max val to be plotted
        self.ymin = float(yrange[0])  # y-axis min val to be plotted
        self.ymax = float(yrange[1])  # y-axis max val to be plotted

        if args.xmin:
            self.xmin = float(args.xmin)
        if args.xmax:
            self.xmax = float(args.xmax)
        if args.ymin:
            self.ymin = float(args.ymin)
        if args.ymax:
            self.ymax = float(args.ymax)

        # Acceptable score limits
        min_val, max_val, min_range = -10, 10, 0.01
        self.xmin, self.xmax = BoundValuesToRange(self.xmin, self.xmax, min_val, max_val, min_range)
        self.ymin, self.ymax = BoundValuesToRange(self.ymin, self.ymax, min_val, max_val, min_range)

        self.plot_title = ""          # graph title
        if args.title:
            self.plot_title = args.title

        self.dot_size = 1             # dot size of scatterplot
        if args.dot_size:
            self.dot_size = float(args.dot_size)

        
    def GetDefaultAxisRanges(self):
        axis_range = {}
        
        axis_range['GD1'] = [ 1.0, 1.8]
        axis_range['GD2'] = [ 1.0, 1.8]
        axis_range['GD3'] = [-0.2, 0.6]
        axis_range['AF1'] = [-1.5, 2.5]
        axis_range['AF2'] = [-1.5, 2.5]
        axis_range['AF3'] = [-1.0, 3.5]
        axis_range['EU1'] = [-1.5, 3.5]
        axis_range['EU2'] = [-1.5, 2.5]
        axis_range['EU3'] = [-2.0, 2.0]
        axis_range['EA1'] = [-1.5, 2.5]
        axis_range['EA2'] = [-2.0, 2.0]
        axis_range['EA3'] = [-2.0, 2.0]
        axis_range['EA4'] = [-1.5, 2.5]
        axis_range['SA1'] = [-1.5, 2.5]
        axis_range['SA2'] = [-1.5, 2.5]
        axis_range['IC1'] = [-2.5, 2.0]
        axis_range['IC2'] = [-1.5, 4.0]
        axis_range['IC3'] = [-2.5, 3.5]
        
        return axis_range        

    def ValidScores(self, args):
        all_scores = ['GD1', 'GD2', 'GD3', 
                      'AF1', 'AF2', 'AF3',
                      'EU1', 'EU2', 'EU3',
                      'EA1', 'EA2', 'EA3', 'EA4',
                      'SA1', 'SA2', 
                      'IC1', 'IC2', 'IC3']

        xvalid, yvalid = False, False
        
        if args.x_score in all_scores:
            xvalid = True

        if args.y_score in all_scores:
            yvalid = True

        error = ""            
        if not xvalid or not yvalid:
            error = "ERROR:"
            if not xvalid:
                error += "\nInvalid x_score: " + str(args.x_score)

            if not yvalid:
                error += "\ninvalid y_score: " + str(args.y_score)

            error = error + "\nValid scores are "
            for score in all_scores:
                error = error + score
                if score != "IC3":
                    error += ", "
            
        return error

def BoundValuesToRange(v1, v2, minv, maxv, min_range):
    """
     Given v1, v2, bound them to [minv, maxv], and set v2 - v1 to a min_range 
    """
    newv1 = v1 if v1 > minv else minv
    newv2 = v2 if v2 > minv else minv

    if newv1 > maxv:
        newv1 = maxv 

    if newv2 > maxv:
        newv2 = maxv 

    if newv2 - newv1 < min_range:
        newv2 = newv1 + min_range

    return newv1, newv2
  
