import os

if __name__ == '__main__':
    # check if param_list.txt exists
    if os.path.exists('param_list.txt'):
        os.remove('param_list.txt')
    
    f = open('param_list.txt', 'w')
    
    p_vals = [0.001]
    lower_bounds = [0.1,0.15,0.2,0.25,0.3]
    upper_bounds = [0.101,0.1501,0.201,0.2501,0.301,0.5,0.75,0.9]
    
    for p in p_vals:
        for l in lower_bounds:
            for u in upper_bounds:
                if l >= u:
                    continue
                f.write(str(l) + ' ' +str(u) + ' ' + str(p) + '\n')