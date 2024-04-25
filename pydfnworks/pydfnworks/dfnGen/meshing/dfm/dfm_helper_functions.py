

def dfm_get_box_domain(domain):

    box_domain = {"x0": None, "x0": None,
                  "y0": None, "y1": None, 
                  "z0": None, "z1": None 
                  }

    # Extent of domain
    box_domain['x0'] = - 0.5*domain['x']
    box_domain['x1'] = 0.5*domain['x'] 
    box_domain['y0'] = - 0.5*domain['y'] 
    box_domain['y1'] = 0.5*domain['y'] 
    box_domain['z0'] = - 0.5*domain['z'] 
    box_domain['z1'] = 0.5*domain['z'] 

    return box_domain 

