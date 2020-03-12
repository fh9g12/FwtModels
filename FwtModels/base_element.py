class BaseElement:
    def CalcKE(self,p):
        raise NotImplementedError("CalcKE not implemented in the current element")
    def CalcPE(self,p):
        raise NotImplementedError("CalcPE not implemented in the current element")
        
