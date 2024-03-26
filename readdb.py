# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 05:26:25 2024

@author: schlammi
"""
import sqlite3

connection = sqlite3.connect('k2viewer.db')
cursor = connection.cursor()
rows = cursor.execute("SELECT run,value,uncertainty, title FROM k2data;").fetchall()

#rows = cursor.execute("""SELECT run,value,uncertainty, 
#                      title FROM k2data WHERE run="{0}";""".format('240325D')).fetchall()
print(rows)
connection.close()
            