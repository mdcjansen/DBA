#!/usr/bin/env python3

from datetime import datetime


date = datetime.now()
#d = date.strftime("%d/%m/%Y %H:%M:%S")
print("date {d}".format(d=date.strftime("%d/%m/%Y %H:%M:%S")))

# print(datetime.strftime("%d/%m/%Y %H:%M:%S"))