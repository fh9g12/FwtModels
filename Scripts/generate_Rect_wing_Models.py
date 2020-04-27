import sys, os
sys.path.insert(1, os.path.join(sys.path[0], '..'))
from FwtModels.RectWing.Models import GenRectWingModel, AeroModel
import multiprocessing as mp

def Model_0(b_modes,t_modes):
    return GenRectWingModel(b_modes,t_modes,fwt_free=False,iwt=True,iwb=False,aero_model=AeroModel.LiftOnly)
def Model_1(b_modes,t_modes):
    return GenRectWingModel(b_modes,t_modes,fwt_free=True,iwt=False,iwb=False,aero_model=AeroModel.LiftOnly)
def Model_2(b_modes,t_modes):
    return GenRectWingModel(b_modes,t_modes,fwt_free=True,iwt=True,iwb=False,aero_model=AeroModel.LiftOnly)
def Model_3(b_modes,t_modes):
    return GenRectWingModel(b_modes,t_modes,fwt_free=True,iwt=True,iwb=True,aero_model=AeroModel.LiftOnly)
def Model_4(b_modes,t_modes):
    return GenRectWingModel(b_modes,t_modes,fwt_free=True,iwt=True,iwb=True,aero_model=AeroModel.LiftOnly_rot)
def Model_5(b_modes,t_modes):
    return GenRectWingModel(b_modes,t_modes,fwt_free=True,iwt=True,iwb=True,aero_model=AeroModel.LiftAndDrag_rot)
def Model_6(b_modes,t_modes):
    model,p = Model_4(b_modes,t_modes)
    return model.linearise(p).msubs({i:0 for i in p.fp}),p
def Model_7(b_modes,t_modes):
    model,p = Model_5(b_modes,t_modes)
    return model.linearise(p).msubs({i:0 for i in p.fp}),p
def Model_8(b_modes,t_modes):
    model,p = GenRectWingModel(b_modes,t_modes,fwt_free=True,iwt=True,iwb=False,aero_model=AeroModel.LiftOnly_rot)
    return model.linearise(p).msubs({i:0 for i in p.fp}),p

ModelFactory = {
    0:Model_0,
    1:Model_1,
    2:Model_2,
    3:Model_3,
    4:Model_4,
    5:Model_5,
    6:Model_6,
    7:Model_7,
    8:Model_8
}

def MakeModel(k,v,b_modes,t_modes):
    print(f'Generating Model {k}')
    try:
        model,_ = v(b_modes,t_modes)
        model.to_file(f'{b_modes}B{t_modes}T-M{k}.py')
        print(f'Generated Model {k}')
    except:
        print(f'Model {k} exited with error')
    return k


b_modes = 3
t_modes = 3

pool = mp.Pool(mp.cpu_count())

for k,v in ModelFactory.items():
    pool.apply_async(MakeModel,args=(k,v,b_modes,t_modes))
pool.close()
pool.join()