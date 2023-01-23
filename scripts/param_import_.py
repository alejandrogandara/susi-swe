#getVariablesFromFile
from pandas import read_excel
import traceback
from datetime import datetime

#import operator
#import json
#from functools import reduce  # forward compatibility for Python 3
# p is the parameter file

#get variables
def getVariablesFromFile(file, site, var=False):
    # var [True, False]
    # file ['file path']
    vars = {}
    try: p = read_excel(file)
    except: 
        print('file does not exists')
        return
    if var:
        try: p=p[(p['commented'] != '#') & (p['path'].isna()) & (p['var'] == var)]
        except: 
            print('variable not found')
            return
    else: p=p[(p['commented'] != '#') & (p['path'].isna())]
    p.reset_index()
    for ind in p.index:
        if p['type'][ind] == 'date': 
            
            date = datetime.strptime(p[site][ind], '%Y-%m-%d %H:%M:%S')
            globals()[p['var'][ind]] = date
            vars[p['var'][ind]] = date
            #datetime.datetime(date.year,date.month,date.day)
        else: 
            globals()[p['var'][ind]] = p[site][ind] 
            vars[p['var'][ind]] = p[site][ind]
        print('variable updated: ' + p['var'][ind]+' = '+ str(p[site][ind]))
    return vars


def getDicInFile(file, site, var=False):
    try: p = read_excel(file)
    except: 
        print('file does not exists')
        return
    #p = p.where(notnull(p), '')
    if var:
        try: p=p[(p['commented'] != '#') & (p['path'].notnull()) & (p['var'] == var)]
        except: 
            print('variable not found')
            return
    else: p=p[(p['commented'] != '#') & (p['path'].notnull())]
    p.reset_index()

    for uvar in p['var'].unique():
        sp = p[p['var'] == uvar]
        print('------------------------------')
        print('dictionary: '+ uvar + '  keys:' + str(len(sp)))
        globals()[uvar]
        for ind in sp.index:
            try: 
                path = p['path'][ind].replace(', ','][')
                if p['type'][ind] == 'str': 
                    value = "'"+ str(p[site][ind]) + "'"
                else: value = str(p[site][ind])
                #exec(p['var'][ind] + p['path'][ind].replace(',','][').replace(' ', '') + ' = ' + str(p[site][ind]))
                line = uvar + path + ' = ' + value
                print(line)
                exec(line)
            #except: print('ERROR: '+ uvar + path + ' = ' + value + ' NOT UPDATED' )
            except Exception: traceback.print_exc()