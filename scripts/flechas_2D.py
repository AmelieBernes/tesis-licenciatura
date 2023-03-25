#Script usado para dibujar en el axis de una figura los ejes x y y con flechas en las ramas positivas.

#Gracias a Felix Hoffman. https://stackoverflow.com/questions/17646247/how-to-make-fuller-axis-arrows-with-matplotlib

def dibujar_flechas_2d(fig, axis):
    xmin, xmax = axis.get_xlim()
    ymin, ymax = axis.get_ylim()
    # get width and height of axes object to compute
    # matching arrowhead length and width
    dps = fig.dpi_scale_trans.inverted()
    bbox = axis.get_window_extent().transformed(dps)
    width, height = bbox.width, bbox.height
    
    # manual arrowhead width and length
    hw = 1./50.*(ymax-ymin)
    hl = 1./50.*(xmax-xmin)
    lw = 1. # axis line width
    ohg = 0.3 # arrow overhang
    
    # compute matching arrowhead length and width
    yhw = hw/(ymax-ymin)*(xmax-xmin)* height/width
    yhl = hl/(xmax-xmin)*(ymax-ymin)* width/height
    
    # draw x and y axis
    axis.arrow(xmin, 0, xmax-xmin, 0., fc='k', ec='k', lw = lw,
             head_width=hw, head_length=hl, overhang = ohg,
             length_includes_head= True, clip_on = False)
    
    axis.arrow(0, ymin, 0., ymax-ymin, fc='k', ec='k', lw = lw,
             head_width=yhw, head_length=yhl, overhang = ohg,
             length_includes_head= True, clip_on = False)
