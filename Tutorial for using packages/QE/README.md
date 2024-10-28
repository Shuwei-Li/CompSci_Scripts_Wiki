There are some tips for learning and using QE:
=============================================

    1. The detailed introduction of parameters used in QE can be found in its official Website: 
    ------------------------------------------------------------------------------------------
    
        https://www.quantum-espresso.org/documentation/input-data-description/
    
    2. A great tutorial video can be found on YouTube: 
    --------------------------------------------------
    
        https://www.youtube.com/watch?v=NKumlrnfrWU&list=PLGntAYRT8AVmQMyurFoncyOdHljqeGU_R&index=3
    
    and the corresponding files can be found on GitHub: 
    --------------------------------------------------
    
        https://github.com/quantumNerd/Quantum-Espresso-Tutorial-2019-Projects
        
    3. The difference between a calculation='bands' and calculation='nscf'：
    -----------------------------------------------------------------------
    
        （1）the former uses exclusively the k-point provided while the latter might add points to respect crystal symmetry;
        （2）In some cases, pw.x calculates additional k-points which are not provided in the k-point list of the input. if this happens, you need to use the keyword of calculation='bands' instead of calculation='nscf'
