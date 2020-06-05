import pandas as pd
import numpy as np

def iterable(obj):
    try:
        iter(obj)
    except Exception:
        return False
    else:
        return True
    
def SeriesFilter(series,var):
    if var is None:
        return series == series
    elif iterable(var):
        return series.isin(var)
    else:
        return series == var

def GetReducedDataSet(data,filter_dict):
    d = data.copy()
    for k,v in filter_dict.items():
        d = d[SeriesFilter(d[k],v)]
    return d

def PlotLines(data,x,y,ax,hue=None,style=None,sigma=5,legend=True):
    styles = ['-','--','-.',':']  
    for style_index,i in enumerate(data[style].unique() if style is not None else [0]):
        ax.set_prop_cycle(None)
        data_style = data[data[style]==i] if style is not None else data
        for hue_index,j in enumerate(data_style[hue].unique()if hue is not None else [0]):
            data_test = data_style[data_style[hue]==j] if hue is not None else data_style
            yi = data_test[y].to_numpy()            
            #clean up jumps in data
            ind = []
            for k in range(0,len(yi)):
                delta_1 = np.abs(yi[k]-yi[k-1])
                delta_2 = np.abs(yi[k-2]-yi[k-1])       
                if (delta_2*sigma <delta_1):
                    ind.append(k-1)
            for index in ind:
                yi[index] = np.NaN

            v = data_test[x].to_numpy()
            for k in range(0,len(v)):
                if v[k]<v[k-1]:
                    v[k-1] = np.NaN 
                elif v[k]-v[k-1]>2:
                    v[k-1] = np.NaN 
            ax.plot(v,yi,styles[style_index])
    if legend:
        ax.legend(data[hue].unique(),title=hue)
        
def FlutterPlot(data,x,y,Modes,hue,ax,style='-',sigma=5,legend=True):
    legend = 'full' if legend else legend
    for i in data[Modes].unique():
        dat = data[data[Modes]==i]
        PlotLines(dat,x,y,ax,hue,style,sigma,legend)
        legend= False
            
            
def CoastDeltaPlot(data,hue, ax,Delta = True):
    # plot Coast Angle Changes
    dataMode0 =data[data['Mode']==0]
    unique_vars = np.unique(data[hue])
    
    coast_data = []
    
    for i in unique_vars:
        a = data[data[hue]==i]['Coast Angle [Deg]'].to_numpy()
        v = data[data[hue]==i]['V'].to_numpy()
        for i in range(1,len(v)):
            if v[i]<v[i-1]:
                v[i-1] = np.NaN
        coast_data.append((v,a))
    
    for v,a in coast_data:
        ax.plot(v,a)
    ax.set_ylabel('Coast Angle [Deg]')
    ax.set_xlabel('Velocity [m/s]')
    if Delta:
        vs = np.linspace(0,max([np.max(v) for v,a in coast_data]),100)   #np.max(v_1),100)
        a_interp = [np.interp(vs,v,a) for v,a in coast_data]
        delta = a_interp[0] - a_interp[1]
        ax2 = ax.twinx()
        ax2.set_ylabel('Delta Coast Angle [Deg]',color='r')
        ax2.plot(vs,delta,'r--')
        ax2.tick_params(axis='y', labelcolor='r')

#def DeltaPlot(data,x,y,hue, ax):
#    # plot Coast Angle Changes
#    dataMode0 =data[data['Mode']==0]
#    unique_vars = np.unique(data[hue])
    
#    coast_data = []
    
#    for i in unique_vars:
#        a = data[data[hue]==i]['Coast Angle [Deg]'].to_numpy()
#        v = data[data[hue]==i]['V'].to_numpy()
#        for i in range(1,len(v)):
#            if v[i]<v[i-1]:
#                v[i-1] = np.NaN
#        coast_data.append((v,a))
    
#    for v,a in coast_data:
#        ax.plot(v,a)
#    ax.set_ylabel('Coast Angle [Deg]')
#    ax.set_xlabel('Velocity [m/s]')
#    if Delta:
#        vs = np.linspace(0,max([np.max(v) for v,a in coast_data]),100)   #np.max(v_1),100)
#        a_interp = [np.interp(vs,v,a) for v,a in coast_data]
#        delta = a_interp[0] - a_interp[1]
#        ax2 = ax.twinx()
#        ax2.set_ylabel('Delta Coast Angle [Deg]',color='r')
#        ax2.plot(vs,delta,'r--')
#        ax2.tick_params(axis='y', labelcolor='r')
