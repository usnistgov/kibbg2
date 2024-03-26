# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 05:26:25 2024

@author: schlammi
"""

import sqlite3
connection = sqlite3.connect('k2viewer.db')
cursor = connection.cursor()
cursor.execute("""
               CREATE TABLE k2data (run TEXT UNIQUE, value FLOAT,
            uncertainty FLOAT, title TEXT, time TIMESTAMP, airdens FLOAT)""")
               
connection.close()
            