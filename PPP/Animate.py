#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Joseph Elmes : NERC-Funded PhD Researcher in Applied Mathematics
University of Leeds : Leeds LS2 9JT : ml14je@leeds.ac.uk

Python 3.7 : Tue Nov  5 15:07:23 2019
"""

import numpy as np
import os

def animate_images(animation_title, total_frames, fps_rate, folder_name):
    from PIL import Image
    import cv2

    for i in range(total_frames):
        file_path = os.path.join(folder_name, f'{i}.png')
        im = np.array(Image.open(file_path))
        h, w = im.shape[0], im.shape[1]
        #print(h,w)
        
        if i == 0:
            h_main, w_main = h, w
            fourcc = cv2.VideoWriter_fourcc(*'MP4V')
            video = cv2.VideoWriter(animation_title+'.mp4', fourcc, fps_rate,
                                    (w_main, h_main))
            
        else:
            assert h == h_main
            assert w == w_main
        video.write(np.array(cv2.cvtColor(im, cv2.COLOR_BGR2RGB)))
    video.release()
    cv2.destroyAllWindows()
    

if __name__ == '__main__':
    folder_name = 'single dose 3mg 90and10'

    path, dirs, files = next(os.walk(folder_name))
    total_frames = len(files)
    animate_images('Cell Evolution', total_frames, 0.1, folder_name)
    
