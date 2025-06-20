This repo contains a MATLAB class containing forward and inverse kinematics for a Comau manipulator.   
The class has animation methods, which run very slow on MATLAB R2025a, due to the new WebGL based graphics system.   

Forward kinematics is generated using screws and their exponential mapping:
$$\mathbf{T}(\mathbf{q})=\mathbf{e}^{\tilde{\mathcal{S}}_1q_1}\mathbf{e}^{\tilde{\mathcal{S}}_2q_2}\mathbf{e}^{\tilde{\mathcal{S}}_3(q_2+q_3)}\mathbf{T}(\mathbf{0})$$   
where    
$$\mathbf{q}\in\mathbb{R}^3,\hspace{1cm}\mathcal{S}=\begin{bmatrix}\mathbf{u} \\ \mathbf{v}\end{bmatrix}\in\mathbb{R}^6,\hspace{1cm}\tilde{\mathcal{S}} = \begin{bmatrix}\tilde{\mathbf{u}} & \mathbf{v}\\\mathbf{0}&0\end{bmatrix}\in se(3)\sub\mathbb{R}^{4\times 4}$$   
where   
$$\tilde{\mathbf{u}}=\begin{bmatrix}0&-u_3&u_2\\u_3&0&-u_1\\-u_2&u_1&0\end{bmatrix}\in so(3)$$   
The vector $\mathbf{u}$ is a unit vector along the screw axis and the vector $\mathbf{p}$ is defined   
$$\mathbf{p}=-\mathbf{u}\times\mathbf{r}$$   
where $\mathbf{r}$ is a vector from the origin to any point along the screw axis.