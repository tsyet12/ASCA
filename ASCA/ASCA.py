import numpy as np
from sklearn.base import TransformerMixin, RegressorMixin, BaseEstimator
import matplotlib.pyplot as plt
from matplotlib.markers import MarkerStyle
import matplotlib.cm as cm
from scipy.spatial import ConvexHull
import math
class ASCA(BaseEstimator):
    def __init__(self):
        self.__name__='ASCA'
        
        #Default
        self.factors=None
        self.interactions=None
        self.data=None
        self.design=None
        self.effects=None
        
        #Total (average)
        self.residuals=None
        self.total_factors=None
        self._total_interactions=None
        
    @staticmethod
    def svd_signstable(X):
        ###Based on Bro, R., Acar, E. and Kolda, T.G., 2008. Resolving the sign ambiguity in the singular value decomposition. Journal of Chemometrics: A Journal of the Chemometrics Society, 22(2), pp.135-140.
        try:
            X=np.asarray(X)
        except:
            pass    
        U, D, V= np.linalg.svd(X,full_matrices=False)
        
        V=V.T  #python V is transposed compared to matlab
        K=len(D)
        s_left=np.zeros((1,K))
        
        #step 1
        for k in range(K):
            select=np.setdiff1d(list(range(K)),k)
            DD=np.zeros((K-1,K-1))
            np.fill_diagonal(DD,D[select])
            Y=X-U[:,select]@DD@V[:,select].T
            
            s_left_parts= np.zeros((1,Y.shape[1]))
            
            for j in range(Y.shape[1]):
                temp_prod=(U[:,k].T)@(Y[:,j])
                s_left_parts[:,j]=(np.sign(temp_prod)+ (temp_prod == 0))*(temp_prod**2)
            
            s_left[:,k]=np.sum(s_left_parts)
           
        #step 2
        s_right=np.zeros((1,K)) 
        for k in range(K):
            select=np.setdiff1d(list(range(K)),k)
            DD=np.zeros((K-1,K-1))
            np.fill_diagonal(DD,D[select])
            Y=X-U[:,select]@DD@V[:,select].T
            
            s_right_parts=np.zeros((1,Y.shape[0]))
            for i in range(Y.shape[0]):
                temp_prod= (V[:,k].T)@(Y[i,:].T)
                s_right_parts[:,i]=(np.sign(temp_prod)+(temp_prod==0))*(temp_prod**2)
            s_right[:,k]=np.sum(s_right_parts)    
         
        #step 3
        for k in range(K):
            if (s_right[:,k]*s_left[:,k])<0:
                if s_left[:,k]<s_right[:,k]:
                    s_left[:,k]=-s_left[:,k]
                else:
                    s_right[:,k]=-s_right[:,k]
        left=np.zeros((K,K))
        right=np.zeros((K,K))
        np.fill_diagonal(left,np.sign(s_left)+ (s_left == 0))
        np.fill_diagonal(right,np.sign(s_right)+ (s_right == 0))
        U=U@left
        V=V@right
        return U, D, V
    
    @classmethod
    def do_PCA(self,_X,_Res):
        u,sv,v=self.svd_signstable(_X)
        #u,sv,v=np.linalg.svd(_X,full_matrices=False)
        scores=u*sv
        singular_values=np.zeros((len(sv), len(sv)))
        np.fill_diagonal(singular_values, sv)
        #singular_values=np.diag(sv)
        loadings = v 
        projected=_Res@v
        explained=(sv**2)/np.sum(sv**2)*100
        return scores, loadings, projected, singular_values, explained
        
    
    def fit(self, X, y, interactions=[None]):
        #initialize for SCA/PCA
        #scores
        self.factors_scores=[]
        self.interaction_scores=[]
        #loadings
        self.factors_loadings=[]
        self.interaction_loadings=[]
        #projected
        self.factors_projected=[]
        self.interaction_projected=[]
        #singular
        self.factors_singular=[]
        self.interaction_singular=[]
        #explained (variance)
        self.factors_explained=[]
        self.interaction_explained=[]
     
        #Save input data
        self.data=X
        self.design=y
        
        #Mean center data
        Xmean=np.mean(X,axis=0)
        Xm=X-Xmean
        #Save as F (sklearn syntax issue)
        F=y
        
        #Make a set of unique factors
        factor_set=set(F.flatten())
        
        #Prepare a zero matrix corresponding to Xm
        zero=np.zeros_like(Xm)
        
        #Calculate Effects
        X_effect=[]        
        for effect in range(F.shape[1]):
          X_effect.append(np.zeros_like(Xm))
          for f in factor_set:
              select=F[:,effect]==f
              select_mat=np.where([select,select],Xm.T,zero.T).T
              select_avg=np.sum(select_mat,axis=0)/(sum(select)+np.finfo(float).eps) #average avoiding empty entry
              select_mat[select_mat.nonzero()]=1 #set all as 1
              X_effect[effect]=X_effect[effect]+select_mat*select_avg
        
        #Sum the effects
        Total_effect=[]
        for Xe in X_effect:
            if len(Total_effect)==0:
                Total_effect=Xe
            else:
                Total_effect=Total_effect+Xe    
        
        #Calculate Interactions
        X_interact=[]
        for inter in range(len(interactions)):
          X_interact.append(np.zeros_like(Xm))
          for f1 in factor_set:
            for f2 in factor_set:
              select=np.sum(F==[f1,f2],axis=1)==2
              select_mat=np.where([select,select],Xm.T,zero.T).T
              select_avg=np.sum(select_mat,axis=0)/(sum(select)+np.finfo(float).eps) #average avoiding empty entry
              select_mat[select_mat.nonzero()]=1 #set all as 1
              X_interact[inter]=X_interact[inter]+select_mat*select_avg-select_mat*Total_effect
              
        #Sum the Interactions
        Total_interact=[]
        for Xi in X_interact:
            if len(Total_interact)==0:
                Total_interact=Xi
            else:
                Total_interact=Total_interact+Xi

        #Calculate Residual/Error
        E=Xm-np.mean(Xm,axis=0)-Total_effect-Total_interact
        
        #Percentage Effects calculation
        SSQ_X=np.sum(X*X)
        SSQ_mean=np.sum(np.tile(Xmean,(X.shape[0],1))**2)
        
        SSQ_factors=[np.sum(Xe**2) for Xe in X_effect] #for each individual factors

        SSQ_residuals=np.sum(E**2)
        SSQ_interactions=np.sum(Total_interact**2)
        SSQ=np.asarray([SSQ_mean,*SSQ_factors,SSQ_interactions,SSQ_residuals])
        percentage_effect=SSQ/SSQ_X*100

        #Save all information
        self.factors=X_effect
        self.interactions=X_interact
        self.total_factors=Total_effect
        self._total_interactions=Total_interact
        self.effects=percentage_effect
        self.residuals=E
        
        #SCA
        for Xe in X_effect:
            _scores, _loadings, _projected, _singular_values, _explained = self.do_PCA(Xe,E)
            self.factors_scores.append(_scores)
            self.factors_loadings.append(_loadings)
            self.factors_projected.append(_projected)
            self.factors_singular.append(_singular_values)
            self.factors_explained.append(_explained)
        
        for Xi in X_interact:
            _scores, _loadings, _projected, _singular_values, _explained = self.do_PCA(Xi,E)
            self.interaction_scores.append(_scores)
            self.interaction_loadings.append(_loadings)
            self.interaction_projected.append(_projected)
            self.interaction_singular.append(_singular_values)
            self.interaction_explained.append(_explained)
        self.factors_scores[1][:,1]=-self.factors_scores[1][:,1] 
    def plot_interactions(self):
        for ii in range(len(self.interaction_scores)):
            #scores
            fig, ax= plt.subplots(1,2)
            ax[0].scatter(self.interaction_scores[ii][:,0],self.interaction_scores[ii][:,1],color='r')
            ax[0].scatter(self.interaction_projected[ii][:,0],self.interaction_projected[ii][:,1],c='blue', marker=MarkerStyle('o',fillstyle='none'))
            n=[" ".join(list(map(str,x))) for x in self.design]

            for i, txt in enumerate(n):
                ax[0].annotate(txt, (self.interaction_scores[ii][i,0], self.interaction_scores[ii][i,1]))
            ax[0].set_xlabel('PC 1 ('+str(round(self.interaction_explained[ii][0],2))+ '% expl. variance)')    
            ax[0].set_ylabel('PC 2 ('+str(round(self.interaction_explained[ii][1],2))+ '% expl. variance)') 
            ax[0].set_title('Interaction='+str(ii+1))
            
            #loadings
            first_loading=self.interaction_loadings[ii][:,0]
            second_loading=self.interaction_loadings[ii][:,1]
            labels=['Response '+str(_i+1) for _i in range(self.interaction_loadings[ii].shape[0])]
            x = np.arange(len(labels))
            width = 0.2
            my_cmap=list(plt.get_cmap("Set1").colors)
            rects1 = ax[1].bar(x - width/2-0.01, first_loading, width, label='First Loading',color=my_cmap[0])
            rects2 = ax[1].bar(x + width/2+0.01, second_loading, width, label='Second Loading',color=my_cmap[1])
            ax[1].hlines(y=0, xmin=0- width-0.05, xmax=len(labels)-1 + width+0.01*5, linewidth=0.5,linestyles='--', color='black')
            ax[1].set_ylabel('Value')
            ax[1].set_title('Loading for Factor '+ str(ii+1))
            ax[1].set_xticks(x, labels)
            ax[1].legend()
            plt.tight_layout()
            plt.show()
        
    def plot_factors(self):
        for ii in range(len(self.factors_scores)):
            #scores
            fig, ax= plt.subplots(1,2)
            
            yy=self.design[:,ii]
            set_yy=list(set(yy))
            colors=list(plt.get_cmap("Set1").colors)
            for _yy in set_yy:
                plot_x=[]
                plot_y=[]
                plot_colors=[]
                plot_label=[]
                
                plot_sx=[]
                plot_sy=[]
                for ss in range(self.factors_projected[ii].shape[0]):
                    if yy[ss]==_yy:
                        cindex=set_yy.index(yy[ss])
                        plot_x.append(self.factors_scores[ii][ss,0]+self.factors_projected[ii][ss,0])
                        plot_y.append(self.factors_scores[ii][ss,1]+self.factors_projected[ii][ss,1])
                        plot_colors.append(colors[cindex])
                        

                        plot_sx=self.factors_scores[ii][ss,0]
                        plot_sy=self.factors_scores[ii][ss,1]
                    else:
                        pass
                ax[0].scatter(plot_sx,plot_sy,marker=MarkerStyle('+'),color=np.array(plot_colors[0]),label='Score level = ' +str(_yy))  
                ax[0].scatter(plot_x,plot_y,c=np.array(plot_colors), marker=MarkerStyle('o',fillstyle='none'),label='Level = '+str(_yy))  
                
                #Boundary box
                points=np.stack((plot_x,plot_y),axis=1)
                hull = ConvexHull(points)
                for simplex in hull.simplices:
                    ax[0].plot(points[simplex, 0], points[simplex, 1], color=np.array(plot_colors[0]),linestyle='--',linewidth=0.3)
  
            n=[" ".join(list(map(str,x))) for x in self.design]
            
            '''
            for i, txt in enumerate(n):
                ax[0].annotate(txt, (self.factors_scores[ii][i,0], self.factors_scores[ii][i,1]))
            ''' 
            
            ax[0].set_xlabel('PC 1 ('+str(round(self.factors_explained[ii][0],2))+ '% expl. variance)')    
            ax[0].set_ylabel('PC 2 ('+str(round(self.factors_explained[ii][1],2))+ '% expl. variance)') 
            ax[0].set_title('Factor '+str(ii+1))
            ax[0].legend()
            
            
            #loadings
            first_loading=self.factors_loadings[ii][:,0]
            second_loading=self.factors_loadings[ii][:,1]
            labels=['Response '+str(_i+1) for _i in range(self.factors_loadings[ii].shape[0])]
            x = np.arange(len(labels))
            width = 0.2
            my_cmap=list(plt.get_cmap("Set1").colors)
            rects1 = ax[1].bar(x - width/2-0.01, first_loading, width, label='First Loading',color=my_cmap[0])
            rects2 = ax[1].bar(x + width/2+0.01, second_loading, width, label='Second Loading',color=my_cmap[1])
            ax[1].hlines(y=0, xmin=0- width-0.05, xmax=len(labels)-1 + width+0.01*5, linewidth=0.5,linestyles='--', color='black')
            ax[1].set_ylabel('Value')
            ax[1].set_title('Loading for Factor '+ str(ii+1))
            ax[1].set_xticks(x, labels)
            ax[1].legend()
            plt.tight_layout()
            plt.show()   

            
    def plot_factors_biplot(self):    
        labels=None
        sign = lambda x: x and 1 - 2 * (x < 0) 
        for ii in range(len(self.factors_scores)):
            xs = self.factors_scores[ii][:,0]
            ys = self.factors_scores[ii][:,1]
            coeff=self.factors_loadings[ii]
            n = len(coeff)
            scalex = 1/(xs.max() - xs.min())
            scaley = 1/(ys.max() - ys.min())
            plt.scatter(xs * scalex,ys * scaley, c = self.design[:,ii],cmap='Set1')
            for i in range(n):
                plt.arrow(0, 0, coeff[i,0], coeff[i,1],color = 'black',alpha = 0.3,head_starts_at_zero=False,head_width=0.03)
                if labels is None:
                    plt.text(coeff[i,0]+sign(coeff[i,0])*0.05, coeff[i,1]+sign(coeff[i,1])*0.05, "Var"+str(i+1), color = 'darkgreen', ha = 'center', va = 'center')
                else:
                    plt.text(coeff[i,0]+sign(coeff[i,0])*0.05, coeff[i,1] +sign(coeff[i,1])*0.05, labels[i], color = 'darkgreen', ha = 'center', va = 'center')
            #plt.xlim(-1,1)
            #plt.ylim(-1,1)
            plt.grid(color='grey', linestyle='--', linewidth=0.5)
            plt.xlabel("PC1 ("+str(round(self.factors_explained[ii][0],2))+"%)")
            plt.ylabel("PC2 ("+str(round(self.factors_explained[ii][1],2))+"%)")
            plt.show()
        
        
        
if __name__=='__main__':
    X = [[1.0000,0.6000], 
    [3.0000,0.4000],
    [2.0000,0.7000],
    [1.0000,0.8000],
    [2.0000,0.0100],
    [2.0000,0.8000],
    [4.0000,1.0000],
    [6.0000,2.0000],
    [5.0000,0.9000],
    [5.0000,1.0000],
    [6.0000,2.0000],
    [5.0000,0.7000]]
    X=np.asarray(X)

    F = [[1,     1],
     [1,     1],
     [1,     2],
     [1,     2],
     [1,     3],
     [1,     3],
     [2,     1],
     [2,     1],
     [2,     2],
     [2,     2],
     [2,     3],
     [2,     3]]
    F=np.asarray(F)
    interactions = [[0, 1]]

    ASCA=ASCA()
    ASCA.fit(X,F,interactions)
    ASCA.plot_factors()
    ASCA.plot_interactions()
    ASCA.plot_factors_biplot()
    #print(ASCA.factors_scores)
    #print(ASCA.factors_loadings)
    #print(ASCA.factors_projected)
    #ang=math.radians(90)
    #A=np.asarray([[math.cos(ang),-math.sin(ang)],[math.sin(ang),math.cos(ang)]])
    #print([ f@A for f  in ASCA.factors_loadings])
    '''
    print(ASCA.effects)
    print(ASCA.factors_scores)
    print(ASCA.factors_explained)
    print(ASCA.interaction_scores)
    print(ASCA.interaction_explained)
    '''
    
