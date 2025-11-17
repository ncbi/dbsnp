#!/usr/bin/env python3

import pandas as pd

class PopColor:
    """
    Set population colors
    
    """

    def __init__(self, args):
        self.background_color = 'gainsboro'
        self.color_list = []
        self.GetColorList()

        if args.bg_color:
            self.background_color = args.bg_color

        if args.color_file:
            self.SetColorListFromFile(args.color_file)
    
        if args.color_str:
            self.SetColorListFromString(args.color_str)

    def SetColorListFromString(self, color_str):
        """
        Read colors from a string passed from command line.
        The string should included a list of valid color names separated by comma,
        e.g., red,pink,#FFFF00
        See https://matplotlib.org/stable/gallery/color/named_colors.html
        for valid color names.
        
        """
        self.color_list = []
        colors = color_str.split(",")
        for color in colors:
            self.color_list.append(color)

    def SetColorListFromFile(self, file):
        """
        Read colors from a string passed from a file
        The string should included a list of valid color names separated by comma,
        e.g., red,pink,#FFFF00
        See https://matplotlib.org/stable/gallery/color/named_colors.html
        for valid color names.
        
        """
        self.color_list = []

        df = pd.read_csv(file, header=None)
        numSbjs = len(df)
        
        for i in range(numSbjs):
            color = df[0].iloc[i]
            self.color_list.append(color)

    def GetColorList(self):
        """
        Only color the top 12 populations. Others are colored using background color.
        
        """
        self.color_list = [
            'blue',
            'red',
            'green',
            'orange',
            'purple',
            'cyan',
            'pink',
            'olive',
            'chocolate',
            'tan',
            'cadetblue',
            'gold'
        ]

