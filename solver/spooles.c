#include <spooles/misc.h>
#include <spooles/FrontMtx.h>
#include <spooles/SymbFac.h>
#include <spooles/Iter/Iter.h>
#include <spooles/Graph.h>

#include <string.h>

#define METHODS 10


/*  allInOne.c  */
/*--------------------------------------------------------------------*/
int spooles_lu_solve(double const *values, int const *column_indices,
                     int n_block_indices, int const *block_indices,
                     double const *rhs, double *result)
{
  /*
    ------------------------------------------------------------------
    all-in-one program to solve A X = Y
  
    (1) read in matrix entries for A and form InpMtx object
    (2) read in right hand side for Y entries and form DenseMtx object
    (3) form Graph object, order matrix and form front tree
    (4) get the permutation, permute the front tree, matrix 
        and right hand side and get the symbolic factorization
    (5) initialize the front matrix object to hold the factor matrices
    (6) compute the numeric factorization
    (7) post-process the factor matrices
    (8) compute the solution
    (9) permute the solution into the original ordering
  
    created -- 98jun04, cca
    ------------------------------------------------------------------
  */
  /*--------------------------------------------------------------------*/
  // char            *matrixFileName, *rhsFileName ;
  DenseMtx        *mtxY, *mtxX ;
  Chv             *rootchv ;
  ChvManager      *chvmanager  ;
  SubMtxManager   *mtxmanager  ;
  FrontMtx        *frontmtx ;
  InpMtx          *mtxA ;
  double          droptol = 0.0, tau = 100. ;
  double          cpus[10] ;
  ETree           *frontETree ;
  FILE            *inputFile, *msgFile ;
  Graph           *graph ;
  int             error, ient, irow, jcol, jrhs, jrow, msglvl, ncol, 
                  nedges, nent, neqns, nrhs, nrow, pivotingflag, seed, 
                  symmetryflag, type ;
  int             *newToOld, *oldToNew ;
  int             stats[20] ;
  IV              *newToOldIV, *oldToNewIV ;
  IVL             *adjIVL, *symbfacIVL ;
  
  int             my_i, my_j, my_k;
  /*--------------------------------------------------------------------*/
  /*
    --------------------
    get input parameters
    --------------------
  */
  // if ( argc != 9 ) {
  //   fprintf(stdout, "\n"
  //       "\n usage: %s msglvl msgFile type symmetryflag pivotingflag"
  //       "\n        matrixFileName rhsFileName seed"
  //       "\n    msglvl -- message level"
  //       "\n    msgFile -- message file"
  //       "\n    type    -- type of entries"
  //       "\n      1 (SPOOLES_REAL)    -- real entries"
  //       "\n      2 (SPOOLES_COMPLEX) -- complex entries"
  //       "\n    symmetryflag -- type of matrix"
  //       "\n      0 (SPOOLES_SYMMETRIC)    -- symmetric entries"
  //       "\n      1 (SPOOLES_HERMITIAN)    -- Hermitian entries"
  //       "\n      2 (SPOOLES_NONSYMMETRIC) -- nonsymmetric entries"
  //       "\n    pivotingflag -- type of pivoting"
  //       "\n      0 (SPOOLES_NO_PIVOTING) -- no pivoting used"
  //       "\n      1 (SPOOLES_PIVOTING)    -- pivoting used"
  //       "\n    matrixFileName -- matrix file name, format"
  //       "\n       nrow ncol nent"
  //       "\n       irow jcol entry"
  //       "\n        ..."
  //       "\n        note: indices are zero based"
  //       "\n    rhsFileName -- right hand side file name, format"
  //       "\n       nrow nrhs "
  //       "\n       ..."
  //       "\n       jrow entry(jrow,0) ... entry(jrow,nrhs-1)"
  //       "\n       ..."
  //       "\n    seed -- random number seed, used for ordering"
  //       "\n", argv[0]) ;
  //   return(0) ;
  // }
  
  msglvl = 0 ;
  msgFile = stdout ; 
  
  srand(time(0x0));
        
  type           = SPOOLES_REAL ;
  symmetryflag   = SPOOLES_NONSYMMETRIC ;
  pivotingflag   = SPOOLES_NO_PIVOTING ;
  seed           = rand() ;
  /*--------------------------------------------------------------------*/
  /*
    --------------------------------------------
    STEP 1: read the entries from the input file 
            and create the InpMtx object
    ------------------------------------double *values, int *column_indices,
                       int n_block_indices, int *block_indices,
                       double *result);--------
  */
  // inputFile = fopen(matrixFileName, "r") ;
  // fscanf(inputFile, "%d %d %d", &nrow, &ncol, &nent) ;
  nrow = n_block_indices - 1;
  ncol = n_block_indices - 1;
  nent = 0;
  neqns = nrow ;
  mtxA = InpMtx_new() ;
  InpMtx_init(mtxA, INPMTX_BY_ROWS, type, nent, neqns) ;
  
  for(my_i = 0, my_k = 0; my_i < n_block_indices - 1; ++my_i)
    for(my_j = block_indices[my_i]; my_j < block_indices[my_i + 1]; ++my_j, ++my_k)
      InpMtx_inputRealEntry(mtxA, my_i, column_indices[my_k] - 1, values[my_k]) ;
  
  InpMtx_changeStorageMode(mtxA, INPMTX_BY_VECTORS) ;
  if ( msglvl > 2 ) {
    fprintf(msgFile, "\n\n input matrix") ;
    InpMtx_writeForHumanEye(mtxA, msgFile) ;
    fflush(msgFile) ;
  }
  /*--------------------------------------------------------------------*/
  /*
    -----------------------------------------
    STEP 2: read the right hand side matrix Y
    -----------------------------------------
  */
  // inputFile = fopen(rhsFileName, "r") ;
  // fscanf(inputFile, "%d %d", &nrow, &nrhs) ;
  nrow = n_block_indices - 1;
  nrhs = 1;
  mtxY = DenseMtx_new() ;
  DenseMtx_init(mtxY, type, 0, 0, neqns, nrhs, 1, neqns) ;
  DenseMtx_zero(mtxY) ;
  
  for(my_i = 0; my_i < n_block_indices - 1; ++my_i)
    DenseMtx_setRealEntry(mtxY, my_i, 0, rhs[my_i]);
    
  if ( msglvl > 2 ) {
    fprintf(msgFile, "\n\n rhs matrix in original ordering") ;
    DenseMtx_writeForHumanEye(mtxY, msgFile) ;
    fflush(msgFile) ;
  }
  /*--------------------------------------------------------------------*/
  /*
    -------------------------------------------------
    STEP 3 : find a low-fill ordering
    (1) create the Graph object
    (2) order the graph using multiple minimum degree
    -------------------------------------------------
  */
  graph = Graph_new() ;
  adjIVL = InpMtx_fullAdjacency(mtxA) ;
  nedges = IVL_tsize(adjIVL) ;
  Graph_init2(graph, 0, neqns, 0, nedges, neqns, nedges, adjIVL,
              NULL, NULL) ;
  if ( msglvl > 2 ) {
    fprintf(msgFile, "\n\n graph of the input matrix") ;
    Graph_writeForHumanEye(graph, msgFile) ;
    fflush(msgFile) ;
  }
  frontETree = orderViaMMD(graph, seed, msglvl, msgFile) ;
  if ( msglvl > 2 ) {
    fprintf(msgFile, "\n\n front tree from ordering") ;
    ETree_writeForHumanEye(frontETree, msgFile) ;
    fflush(msgFile) ;
  }
  /*--------------------------------------------------------------------*/
  /*
    ----------------------------------------------------
    STEP 4: get the permutation, permute the front tree,
            permute the matrix and right hand side, and 
            get the symbolic factorization
    ----------------------------------------------------
  */
  oldToNewIV = ETree_oldToNewVtxPerm(frontETree) ;
  oldToNew   = IV_entries(oldToNewIV) ;
  newToOldIV = ETree_newToOldVtxPerm(frontETree) ;
  newToOld   = IV_entries(newToOldIV) ;
  ETree_permuteVertices(frontETree, oldToNewIV) ;
  InpMtx_permute(mtxA, oldToNew, oldToNew) ;
  if (  symmetryflag == SPOOLES_SYMMETRIC
    || symmetryflag == SPOOLES_HERMITIAN ) {
    InpMtx_mapToUpperTriangle(mtxA) ;
  }
  InpMtx_changeCoordType(mtxA, INPMTX_BY_CHEVRONS) ;
  InpMtx_changeStorageMode(mtxA, INPMTX_BY_VECTORS) ;
  DenseMtx_permuteRows(mtxY, oldToNewIV) ;
  symbfacIVL = SymbFac_initFromInpMtx(frontETree, mtxA) ;
  if ( msglvl > 2 ) {
    fprintf(msgFile, "\n\n old-to-new permutation vector") ;
    IV_writeForHumanEye(oldToNewIV, msgFile) ;
    fprintf(msgFile, "\n\n new-to-old permutation vector") ;
    IV_writeForHumanEye(newToOldIV, msgFile) ;
    fprintf(msgFile, "\n\n front tree after permutation") ;
    ETree_writeForHumanEye(frontETree, msgFile) ;
    fprintf(msgFile, "\n\n input matrix after permutation") ;
    InpMtx_writeForHumanEye(mtxA, msgFile) ;
    fprintf(msgFile, "\n\n right hand side matrix after permutation") ;
    DenseMtx_writeForHumanEye(mtxY, msgFile) ;
    fprintf(msgFile, "\n\n symbolic factorization") ;
    IVL_writeForHumanEye(symbfacIVL, msgFile) ;
    fflush(msgFile) ;
  }
  /*--------------------------------------------------------------------*/
  /*
    ------------------------------------------
    STEP 5: initialize the front matrix object
    ------------------------------------------
  */
  frontmtx   = FrontMtx_new() ;
  mtxmanager = SubMtxManager_new() ;
  SubMtxManager_init(mtxmanager, NO_LOCK, 0) ;
  FrontMtx_init(frontmtx, frontETree, symbfacIVL, type, symmetryflag, 
                FRONTMTX_DENSE_FRONTS, pivotingflag, NO_LOCK, 0, NULL, 
                mtxmanager, msglvl, msgFile) ;
  /*--------------------------------------------------------------------*/
  /*
    -----------------------------------------
    STEP 6: compute the numeric factorization
    -----------------------------------------
  */
  chvmanager = ChvManager_new() ;
  ChvManager_init(chvmanager, NO_LOCK, 1) ;
  DVfill(10, cpus, 0.0) ;
  IVfill(20, stats, 0) ;
  rootchv = FrontMtx_factorInpMtx(frontmtx, mtxA, tau, droptol, 
              chvmanager, &error, cpus, stats, msglvl, msgFile) ;
  ChvManager_free(chvmanager) ;
  if ( msglvl > 2 ) {
    fprintf(msgFile, "\n\n factor matrix") ;
    FrontMtx_writeForHumanEye(frontmtx, msgFile) ;
    fflush(msgFile) ;
  }
  if ( rootchv != NULL ) {
    fprintf(msgFile, "\n\n matrix found to be singular\n") ;
    exit(-1) ;
  }
  if ( error >= 0 ) {
    fprintf(msgFile, "\n\n error encountered at front %d", error) ;
    exit(-1) ;
  }
  /*--------------------------------------------------------------------*/
  /*
    --------------------------------------
    STEP 7: post-process the factorization
    --------------------------------------
  */
  FrontMtx_postProcess(frontmtx, msglvl, msgFile) ;
  if ( msglvl > 2 ) {
    fprintf(msgFile, "\n\n factor matrix after post-processing") ;
    FrontMtx_writeForHumanEye(frontmtx, msgFile) ;
    fflush(msgFile) ;
  }
  /*--------------------------------------------------------------------*/
  /*
    -------------------------------
    STEP 8: solve the linear system
    -------------------------------
  */
  mtxX = DenseMtx_new() ;
  DenseMtx_init(mtxX, type, 0, 0, neqns, nrhs, 1, neqns) ;
  DenseMtx_zero(mtxX) ;
  FrontMtx_solve(frontmtx, mtxX, mtxY, mtxmanager, 
                cpus, msglvl, msgFile) ;
  if ( msglvl > 2 ) {
    fprintf(msgFile, "\n\n solution matrix in new ordering") ;
    DenseMtx_writeForHumanEye(mtxX, msgFile) ;
    fflush(msgFile) ;
  }
  
  
  /*--------------------------------------------------------------------*/
  /*
    -------------------------------------------------------
    STEP 9: permute the solution into the original ordering
    -------------------------------------------------------
  */
  DenseMtx_permuteRows(mtxX, newToOldIV) ;
  if ( msglvl > 0 ) {
    fprintf(msgFile, "\n\n solution matrix in original ordering") ;
    DenseMtx_writeForHumanEye(mtxX, msgFile) ;
    fflush(msgFile) ;
  }
  
  for(my_i = 0; my_i < n_block_indices - 1; ++my_i)
    DenseMtx_realEntry(mtxX, my_i, 0, &result[my_i]);
    
  /*--------------------------------------------------------------------*/
  /*
    -----------
    free memory
    -----------
  */
  FrontMtx_free(frontmtx) ;
  DenseMtx_free(mtxX) ;
  DenseMtx_free(mtxY) ;
  IV_free(newToOldIV) ;
  IV_free(oldToNewIV) ;
  InpMtx_free(mtxA) ;
  ETree_free(frontETree) ;
  IVL_free(symbfacIVL) ;
  SubMtxManager_free(mtxmanager) ;
  Graph_free(graph) ;
  /*--------------------------------------------------------------------*/
  return(1) ;
}
/*--------------------------------------------------------------------*/ 




int InpMtx_createGraph(InpMtx *mtxA, Graph    *graph);


/*--------------------------------------------------------------------*/
int
spooles_bicgstab_solve(double const *values, int const *column_indices,
                     int n_block_indices, int const *block_indices,
                     double const *rhs, double *result)
/*
   -----------------------------------------------------
   test the factor method for a grid matrix
   (0) read in matrix from source file 
   (1) conver data matrix to InpMtx object if necessary
   (2) create Graph and ETree object if necessary
   (3) read in/create an ETree object
   (4) create a solution matrix object
   (5) multiply the solution with the matrix
       to get a right hand side matrix object
   (6) factor the matrix 
   (7) solve the system

   created   -- 98dec30, jwu
   -----------------------------------------------------
*/
{
  char            ctemp[81];
  Chv             *chv, *rootchv ;
  ChvManager      *chvmanager ;
  DenseMtx        *mtxB, *mtxQ, *mtxX, *mtxZ ;
  double          one[2] = { 1.0, 0.0 } ;
  FrontMtx        *frontmtx ;
  InpMtx          *mtxA ;
  SubMtxManager   *mtxmanager ;
  double          cputotal, droptol, conv_tol, factorops ;
  double          cpus[9] ;
  Drand           drand ;
  double          nops, tau, t1, t2   ;
  ETree           *frontETree   ;
  Graph           *graph ;
  FILE            *msgFile, *inFile ;
  int             error, loc, msglvl, neqns, nzf, iformat, 
                  pivotingflag, rc, seed, sparsityflag, symmetryflag, 
                  method, type, nrhs, etreeflag ;
  int             stats[6] ;
  int             nnzA, Ik, itermax, zversion, iterout ;
  IV              *newToOldIV, *oldToNewIV ;
  IVL             *symbfacIVL ;
  int             i, j, k, m, n, imethod, maxdomainsize, maxzeros, maxsize;
  int             nouter,ninner ;
    
  int             my_i, my_j, my_k;
  
  etreeflag = 1;
  // msgFile = stdout;
  msgFile = fopen("spooles.log", "w");
  msglvl = 0;
  seed = rand();
  method = 0; // BiCGStabR
  nrhs = 1;
  // Ik
  itermax = 0;
  // iterout
  symmetryflag = SPOOLES_NONSYMMETRIC;
  sparsityflag = FRONTMTX_DENSE_FRONTS;
  pivotingflag = SPOOLES_PIVOTING;
  tau = 1;
  droptol = 1E-3;
  conv_tol = 1E-5;
  type = 1;
  neqns = n_block_indices - 1;
  
    
  /*
    --------------------------------------
    initialize the random number generator
    --------------------------------------
  */
  Drand_setDefaultFields(&drand) ;
  Drand_init(&drand) ;
  Drand_setSeed(&drand, seed) ;
  /*Drand_setUniform(&drand, 0.0, 1.0) ;*/
  Drand_setNormal(&drand, 0.0, 1.0) ;
  /*
    ----------------------------------------------
    read in or convert source to the InpMtx object
    ----------------------------------------------
  */
  rc = 1;

  mtxA = InpMtx_new() ;
  InpMtx_init(mtxA, INPMTX_BY_ROWS, type, 0, neqns) ;
  
  for(my_i = 0, my_k = 0; my_i < n_block_indices - 1; ++my_i)
    for(my_j = block_indices[my_i]; my_j < block_indices[my_i + 1]; ++my_j, ++my_k)
      InpMtx_inputRealEntry(mtxA, my_i, column_indices[my_k] - 1, values[my_k]) ;
        
  /*
    Get the nonzeros in matrix A and print it
    */
  nnzA  = InpMtx_nent( mtxA );
  fprintf(msgFile, "\n\n Input matrix size  %d NNZ  %d",
    neqns, nnzA) ;

  
  /*
    --------------------------------------------------------
    generate the linear system
    1. generate solution matrix and fill with random numbers
    2. generate rhs matrix and fill with zeros
    3. compute matrix-matrix multiply
    --------------------------------------------------------
  */
  MARKTIME(t1) ;
  mtxX = DenseMtx_new() ;
  DenseMtx_init(mtxX, type, 0, -1, neqns, nrhs, 1, neqns) ;
  mtxB = DenseMtx_new() ; 
  
  DenseMtx_init(mtxB, type, 0, -1, neqns, nrhs, 1, neqns) ;
  
  for(my_i = 0; my_i < n_block_indices - 1; ++my_i)
    DenseMtx_setRealEntry(mtxB, my_i, 0, rhs[my_i]);
    
  DenseMtx_zero(mtxX) ;
    
  MARKTIME(t2) ;
  fprintf(msgFile, "\n CPU %8.3f : set up the solution and rhs ",
          t2 - t1) ;
  if ( msglvl > 2 ) {
    fprintf(msgFile, "\n\n original mtxA") ;
    InpMtx_writeForHumanEye(mtxA, msgFile) ;
    fprintf(msgFile, "\n\n original mtxX") ;
    DenseMtx_writeForHumanEye(mtxX, msgFile) ;
    fprintf(msgFile, "\n\n original mtxB") ;
    DenseMtx_writeForHumanEye(mtxB, msgFile) ;
    fflush(msgFile) ;
  }
  if (rc != 1) {
    InpMtx_free(mtxA);
    DenseMtx_free(mtxX);
    DenseMtx_free(mtxB);
    goto end_init;
  }
  
  /*
    ------------------------
    read in/create the ETree object
    ------------------------
  */
  
  MARKTIME(t1) ;

  graph = Graph_new() ;
  rc = InpMtx_createGraph(mtxA, graph);
  if (rc!=1) {
    fprintf(msgFile, "\n return value %d from InpMtx_createGraph(%p,%p)",
      rc, mtxA, graph) ;
    Graph_free(graph);
    goto end_tree;
  }
  if (etreeflag == 1) { /* Via BestOfNDandMS */
    maxdomainsize = 500; maxzeros      = 1000; maxsize       = 64    ;
    frontETree = orderViaBestOfNDandMS(graph, maxdomainsize, maxzeros,
              maxsize, seed, msglvl, msgFile) ;
  }
  else if (etreeflag == 2) { /* Via MMD */
    frontETree = orderViaMMD(graph, seed, msglvl, msgFile) ;        
  }
  else if (etreeflag == 3) { /* Via MS */
    maxdomainsize = 500;
    frontETree = orderViaMS(graph, maxdomainsize, seed, msglvl, msgFile) ;
  }
  else if (etreeflag == 4) { /* Via ND */
    maxdomainsize = 500;
    frontETree = orderViaND(graph, maxdomainsize, seed, msglvl, msgFile) ;
  }
  Graph_free(graph);
  
  MARKTIME(t2) ;
  fprintf(msgFile, "\n CPU %8.3f : read in/create frontETree",
    t2 - t1) ;
  if ( rc != 1 ) {
    ETree_free(frontETree);
    goto end_tree;
  }
  
  ETree_leftJustify(frontETree) ;
  if ( msglvl > 1 ) {
    fprintf(msgFile, "\n\n after creating ETree object") ;
    if ( msglvl == 2 ) {
      ETree_writeStats(frontETree, msgFile) ;
    } else {
      ETree_writeForHumanEye(frontETree, msgFile) ;
    }
  }
  fflush(msgFile) ;
  /*
    --------------------------------------------------
    get the permutations, permute the matrix and the 
    front tree, and compute the symbolic factorization
    --------------------------------------------------
  */
  MARKTIME(t1) ;
  oldToNewIV = ETree_oldToNewVtxPerm(frontETree) ;
  newToOldIV = ETree_newToOldVtxPerm(frontETree) ;
  MARKTIME(t2) ;
  fprintf(msgFile, "\n CPU %8.3f : get permutations", t2 - t1) ;
  MARKTIME(t1) ;
  ETree_permuteVertices(frontETree, oldToNewIV) ;
  MARKTIME(t2) ;
  fprintf(msgFile, "\n CPU %8.3f : permute front tree", t2 - t1) ;
  MARKTIME(t1) ;
  InpMtx_permute(mtxA, IV_entries(oldToNewIV), IV_entries(oldToNewIV)) ;
  MARKTIME(t2) ;
  fprintf(msgFile, "\n CPU %8.3f : permute mtxA", t2 - t1) ;
  if (  symmetryflag == SPOOLES_SYMMETRIC
    || symmetryflag == SPOOLES_HERMITIAN ) {
    MARKTIME(t1) ;
    InpMtx_mapToUpperTriangle(mtxA) ;
    MARKTIME(t2) ;
    fprintf(msgFile, "\n CPU %8.3f : map to upper triangle", t2 - t1) ;
  }
  if ( ! INPMTX_IS_BY_CHEVRONS(mtxA) ) {
    MARKTIME(t1) ;
    InpMtx_changeCoordType(mtxA, INPMTX_BY_CHEVRONS) ;
    MARKTIME(t2) ;
    fprintf(msgFile, "\n CPU %8.3f : change coordinate type", t2 - t1) ;
  }
  if ( INPMTX_IS_RAW_DATA(mtxA) ) {
    MARKTIME(t1) ;
    InpMtx_changeStorageMode(mtxA, INPMTX_SORTED) ;
    MARKTIME(t2) ;
    fprintf(msgFile, "\n CPU %8.3f : sort entries ", t2 - t1) ;
  }
  if ( INPMTX_IS_SORTED(mtxA) ) {
    MARKTIME(t1) ;
    InpMtx_changeStorageMode(mtxA, INPMTX_BY_VECTORS) ;
    MARKTIME(t2) ;
    fprintf(msgFile, "\n CPU %8.3f : convert to vectors ", t2 - t1) ;
  }
  MARKTIME(t1) ;
  symbfacIVL = SymbFac_initFromInpMtx(frontETree, mtxA) ;
  MARKTIME(t2) ;
  fprintf(msgFile, "\n CPU %8.3f : symbolic factorization", t2 - t1) ;
  MARKTIME(t1) ;
  DenseMtx_permuteRows(mtxB, oldToNewIV) ;
  MARKTIME(t2) ;
  fprintf(msgFile, "\n CPU %8.3f : permute rhs", t2 - t1) ;
  
  /*
    ------------------------------
    initialize the FrontMtx object
    ------------------------------
  */
  MARKTIME(t1) ;
  frontmtx   = FrontMtx_new() ;
  mtxmanager = SubMtxManager_new() ;
  SubMtxManager_init(mtxmanager, NO_LOCK, 0) ;
  FrontMtx_init(frontmtx, frontETree, symbfacIVL,
                type, symmetryflag, sparsityflag, pivotingflag,
                NO_LOCK, 0, NULL, mtxmanager, msglvl, msgFile) ;
  MARKTIME(t2) ;
  fprintf(msgFile, "\n\n CPU %8.3f : initialize the front matrix",
          t2 - t1) ;
  if ( msglvl > 1 ) {
    fprintf(msgFile,
            "\n nendD  = %d, nentL = %d, nentU = %d",
            frontmtx->nentD, frontmtx->nentL, frontmtx->nentU) ;
    SubMtxManager_writeForHumanEye(mtxmanager, msgFile) ;
  }
  if ( msglvl > 2 ) {
    fprintf(msgFile, "\n front matrix initialized") ;
    FrontMtx_writeForHumanEye(frontmtx, msgFile) ;
    fflush(msgFile) ;
  }
  /*
    -----------------
    factor the matrix
    -----------------
  */
  nzf       = ETree_nFactorEntries(frontETree, symmetryflag) ;
  factorops = ETree_nFactorOps(frontETree, type, symmetryflag) ;
  fprintf(msgFile, 
          "\n %d factor entries, %.0f factor ops, %8.3f ratio",
          nzf, factorops, factorops/nzf) ;
  IVzero(6, stats) ;
  DVzero(9, cpus) ;
  chvmanager = ChvManager_new() ;
  ChvManager_init(chvmanager, NO_LOCK, 1) ;
  MARKTIME(t1) ;
  rootchv = FrontMtx_factorInpMtx(frontmtx, mtxA, tau, droptol, 
                                  chvmanager, &error, cpus, 
                                  stats, msglvl, msgFile) ;
  MARKTIME(t2) ;
  fprintf(msgFile, "\n\n CPU %8.3f : factor matrix, %8.3f mflops",
          t2 - t1, 1.e-6*factorops/(t2-t1)) ;
  if ( rootchv != NULL ) {
    fprintf(msgFile, "\n\n factorization did not complete") ;
    for ( chv = rootchv ; chv != NULL ; chv = chv->next ) {
        fprintf(stdout, "\n chv %d, nD = %d, nL = %d, nU = %d",
                chv->id, chv->nD, chv->nL, chv->nU) ;
    }
  }
  if ( error >= 0 ) {
    fprintf(msgFile, "\n\n error encountered at front %d\n", error) ;
    rc=error ;
    goto end_front;
  }
  fprintf(msgFile,
          "\n %8d pivots, %8d pivot tests, %8d delayed rows and columns",
          stats[0], stats[1], stats[2]) ;
  if ( frontmtx->rowadjIVL != NULL ) {
    fprintf(msgFile,
            "\n %d entries in rowadjIVL", frontmtx->rowadjIVL->tsize) ;
  }
  if ( frontmtx->coladjIVL != NULL ) {
    fprintf(msgFile,
            ", %d entries in coladjIVL", frontmtx->coladjIVL->tsize) ;
  }
  if ( frontmtx->upperblockIVL != NULL ) {
    fprintf(msgFile,
            "\n %d fronts, %d entries in upperblockIVL", 
            frontmtx->nfront, frontmtx->upperblockIVL->tsize) ;
  }
  if ( frontmtx->lowerblockIVL != NULL ) {
    fprintf(msgFile,
            ", %d entries in lowerblockIVL", 
            frontmtx->lowerblockIVL->tsize) ;
  }
  fprintf(msgFile,
          "\n %d entries in D, %d entries in L, %d entries in U",
          stats[3], stats[4], stats[5]) ;
  fprintf(msgFile, "\n %d locks", frontmtx->nlocks) ;
  if (  FRONTMTX_IS_SYMMETRIC(frontmtx)
    || FRONTMTX_IS_HERMITIAN(frontmtx) ) {
    int   nneg, npos, nzero ;
  
    FrontMtx_inertia(frontmtx, &nneg, &nzero, &npos) ;
    fprintf(msgFile, 
            "\n %d negative, %d zero and %d positive eigenvalues",
            nneg, nzero, npos) ;
    fflush(msgFile) ;
  }
  cputotal = cpus[8] ;
  if ( cputotal > 0.0 ) {
    fprintf(msgFile,
    "\n    initialize fronts       %8.3f %6.2f"
    "\n    load original entries   %8.3f %6.2f"
    "\n    update fronts           %8.3f %6.2f"
    "\n    assemble postponed data %8.3f %6.2f"
    "\n    factor fronts           %8.3f %6.2f"
    "\n    extract postponed data  %8.3f %6.2f"
    "\n    store factor entries    %8.3f %6.2f"
    "\n    miscellaneous           %8.3f %6.2f"
    "\n    total time              %8.3f",
    cpus[0], 100.*cpus[0]/cputotal,
    cpus[1], 100.*cpus[1]/cputotal,
    cpus[2], 100.*cpus[2]/cputotal,
    cpus[3], 100.*cpus[3]/cputotal,
    cpus[4], 100.*cpus[4]/cputotal,
    cpus[5], 100.*cpus[5]/cputotal,
    cpus[6], 100.*cpus[6]/cputotal,
    cpus[7], 100.*cpus[7]/cputotal, cputotal) ;
  }
  if ( msglvl > 1 ) {
    SubMtxManager_writeForHumanEye(mtxmanager, msgFile) ;
    ChvManager_writeForHumanEye(chvmanager, msgFile) ;
  }
  if ( msglvl > 2 ) {
    fprintf(msgFile, "\n\n front factor matrix") ;
    FrontMtx_writeForHumanEye(frontmtx, msgFile) ;
  }
  
  /*
    ------------------------------
    post-process the factor matrix
    ------------------------------
  */
  MARKTIME(t1) ;
  FrontMtx_postProcess(frontmtx, msglvl, msgFile) ;
  MARKTIME(t2) ;
  fprintf(msgFile, "\n\n CPU %8.3f : post-process the matrix", t2 - t1) ;
  if ( msglvl > 2 ) {
    fprintf(msgFile, "\n\n front factor matrix after post-processing") ;
    FrontMtx_writeForHumanEye(frontmtx, msgFile) ;
  }
  fprintf(msgFile, "\n\n after post-processing") ;
  if ( msglvl > 1 ) SubMtxManager_writeForHumanEye(frontmtx->manager, msgFile) ;
  /*
    ----------------
    solve the system
    ----------------
  */
  neqns = mtxB->nrow ;
  mtxZ  = DenseMtx_new() ;
  DenseMtx_init(mtxZ, type, 0, 0, neqns, nrhs, 1, neqns) ;
  zversion=INPMTX_IS_COMPLEX_ENTRIES(mtxA);
  
  DenseMtx_zero(mtxZ) ;
  if ( msglvl > 2 ) {
    fprintf(msgFile, "\n\n rhs") ;
    DenseMtx_writeForHumanEye(mtxB, msgFile) ;
    fflush(stdout) ;
  }
  fprintf(msgFile, "\n\n itemax  %d", itermax) ;
  DVzero(6, cpus) ;
  MARKTIME(t1) ;
  switch ( method ) {
  case BiCGStabR :
    if (zversion)
      rc=zbicgstabr(neqns, type, symmetryflag, mtxA, frontmtx, mtxZ, mtxB,
        itermax, conv_tol, msglvl, msgFile);
    else
      rc=bicgstabr(neqns, type, symmetryflag, mtxA, frontmtx, mtxZ, mtxB,
      itermax, conv_tol, msglvl, msgFile);

    break;
  case BiCGStabL :
    if (zversion)
    rc=zbicgstabl(neqns, type, symmetryflag, mtxA, frontmtx, mtxZ, mtxB,
      itermax, conv_tol, msglvl, msgFile);
    else
      rc=bicgstabl(neqns, type, symmetryflag, mtxA, frontmtx, mtxZ, mtxB,
      itermax, conv_tol, msglvl, msgFile);
    break;
  case TFQMRR :
    if (zversion)
      rc=ztfqmrr(neqns, type, symmetryflag, mtxA, frontmtx, mtxZ, mtxB,
    itermax, conv_tol, msglvl, msgFile);
    else
      rc=tfqmrr(neqns, type, symmetryflag, mtxA, frontmtx, mtxZ, mtxB,
    itermax, conv_tol, msglvl, msgFile);
    break;
  case TFQMRL :
    if (zversion)
      rc=ztfqmrl(neqns, type, symmetryflag, mtxA, frontmtx, mtxZ, mtxB,
    itermax, conv_tol, msglvl, msgFile);
    else
      rc=tfqmrl(neqns, type, symmetryflag, mtxA, frontmtx, mtxZ, mtxB,
    itermax, conv_tol, msglvl, msgFile);
    break;
  case PCGR :
    if (zversion)
      rc=zpcgr(neqns, type, symmetryflag, mtxA, frontmtx, mtxZ, mtxB,
        itermax, conv_tol, msglvl, msgFile);
    else
      rc=pcgr(neqns, type, symmetryflag, mtxA, frontmtx, mtxZ, mtxB,
        itermax, conv_tol, msglvl, msgFile);
    break;
  case PCGL :
    if (zversion)
      rc=zpcgl(neqns, type, symmetryflag, mtxA, frontmtx, mtxZ, mtxB,
        itermax, conv_tol, msglvl, msgFile);
    else
      rc=pcgl(neqns, type, symmetryflag, mtxA, frontmtx, mtxZ, mtxB,
        itermax, conv_tol, msglvl, msgFile);
    break;
  case MLBiCGStabR :
    mtxQ = DenseMtx_new() ;
    DenseMtx_init(mtxQ, type, 0, -1, neqns, Ik, 1, neqns) ;
    Drand_setUniform(&drand, 0.0, 1.0) ;
    DenseMtx_fillRandomEntries(mtxQ, &drand) ;
    if (zversion)
      rc=zmlbicgstabr(neqns, type, symmetryflag, mtxA, frontmtx, mtxQ, mtxZ, 
          mtxB, itermax, conv_tol, msglvl, msgFile);
    else
      rc=mlbicgstabr(neqns, type, symmetryflag, mtxA, frontmtx, mtxQ, mtxZ, 
        mtxB, itermax, conv_tol, msglvl, msgFile);
    DenseMtx_free(mtxQ) ;
    break;
  case MLBiCGStabL :
    mtxQ = DenseMtx_new() ;
    DenseMtx_init(mtxQ, type, 0, -1, neqns, Ik, 1, neqns) ;
    Drand_setUniform(&drand, 0.0, 1.0) ;
    DenseMtx_fillRandomEntries(mtxQ, &drand) ;
    if (zversion)
      rc=zmlbicgstabl(neqns, type, symmetryflag, mtxA, frontmtx, mtxQ, mtxZ, 
          mtxB, itermax, conv_tol, msglvl, msgFile);
    else
      rc=mlbicgstabl(neqns, type, symmetryflag, mtxA, frontmtx, mtxQ, mtxZ, 
        mtxB, itermax, conv_tol, msglvl, msgFile);
    DenseMtx_free(mtxQ) ;
    break;
  case BGMRESR:    
    if (zversion)
      fprintf(msgFile, "\n\n *** BGMRESR complex version is not available "
        "at this moment.   ") ;
    else
      rc=bgmresr(neqns, type, symmetryflag, mtxA, frontmtx, mtxZ,
                mtxB, iterout, itermax, &nouter, &ninner, conv_tol,
                msglvl, msgFile);
    break;
  case BGMRESL:    
    if (zversion)
      fprintf(msgFile, "\n\n *** BGMRESR complex version is not available "
        "at this moment.   ") ;
    else
      rc=bgmresl(neqns, type, symmetryflag, mtxA, frontmtx, mtxZ,
                mtxB, iterout, itermax, &nouter, &ninner, conv_tol,
                msglvl, msgFile);
    break;
  default:
    fprintf(msgFile, "\n\n *** Invalid method number   ") ;
  }
  
  MARKTIME(t2) ;
  fprintf(msgFile, "\n\n CPU %8.3f : solve the system", t2 - t1) ;
  if ( msglvl > 2 ) {
    fprintf(msgFile, "\n\n computed solution") ;
    DenseMtx_writeForHumanEye(mtxZ, msgFile) ;
    fflush(stdout) ;
  }
  
/*
  -------------------------------------------------------------
  permute the computed solution back into the original ordering
  -------------------------------------------------------------
*/
  MARKTIME(t1) ;
  DenseMtx_permuteRows(mtxZ, newToOldIV) ;
  MARKTIME(t2) ;
  fprintf(msgFile, "\n CPU %8.3f : permute solution", t2 - t1) ;
  if ( msglvl > 2 ) {
    fprintf(msgFile, "\n\n permuted solution") ;
    DenseMtx_writeForHumanEye(mtxZ, msgFile) ;
    fflush(stdout) ;
  }
/*
  -------------
  save solution
  -------------
// */
//   if (  strcmp(slnFileName, "none") != 0 ) {
//     DenseMtx_writeToFile(mtxZ, slnFileName) ;
//   }
  for(my_i = 0; my_i < n_block_indices - 1; ++my_i)
    DenseMtx_realEntry(mtxZ, my_i, 0, &result[my_i]);
/*
  -----------------
  compute the error
  -----------------
*/
  DenseMtx_sub(mtxZ, mtxX) ;
  if (method <8) {
    mtxQ = DenseMtx_new() ;
    DenseMtx_init(mtxQ, type, 0, -1, neqns, 1, 1, neqns) ;
    rc=DenseMtx_initAsSubmatrix (mtxQ, mtxZ, 0, neqns-1, 0, 0);
    fprintf(msgFile, "\n\n maxabs error = %12.4e", DenseMtx_maxabs(mtxQ)) ;
    DenseMtx_free(mtxQ) ;
  }
  else
    fprintf(msgFile, "\n\n maxabs error = %12.4e", DenseMtx_maxabs(mtxZ)) ;

  if ( msglvl > 1 ) {
    fprintf(msgFile, "\n\n error") ;
    DenseMtx_writeForHumanEye(mtxZ, msgFile) ;
    fflush(stdout) ;
  }
  if ( msglvl > 1 ) 
    SubMtxManager_writeForHumanEye(frontmtx->manager, msgFile) ;
        
  /*
    ------------------------
    free the working storage
    ------------------------
  */
  DenseMtx_free(mtxZ) ;
  
  end_front:
  ChvManager_free(chvmanager) ;
  SubMtxManager_free(mtxmanager) ;
  FrontMtx_free(frontmtx) ;
  IVL_free(symbfacIVL) ;
  IV_free(oldToNewIV) ;
  IV_free(newToOldIV) ;
  
  end_tree:
  ETree_free(frontETree) ;
  
  end_init:
  DenseMtx_free(mtxB) ;
  DenseMtx_free(mtxX) ;
  
  end_read:
  InpMtx_free(mtxA) ;
  
  fprintf(msgFile, "\n") ;
  fclose(msgFile) ;
  
  return(rc) ;

}


/*--------------------------------------------------------------------*/

/*  createGraph.c  */
/*--------------------------------------------------------------------*/
int
InpMtx_createGraph (InpMtx *mtxA, Graph    *graph)
/*
   ----------------------------------------------------
   read in a InpMtx object and create the Graph object

   created -- 97feb14, cca
   ----------------------------------------------------
*/
{
  int      count, msglvl, nvtx, rc ;
  IVL      *adjIVL ;
  
  InpMtx_changeStorageMode(mtxA, 3) ;
  nvtx  = 1 + IV_max(&mtxA->ivec1IV) ;
  count = 1 + IV_max(&mtxA->ivec2IV) ;
  if ( nvtx < count ) {
    nvtx = count ;
  }
  /*
    ------------------------------------
    create the full adjacency IVL object
    ------------------------------------
  */
  adjIVL = InpMtx_fullAdjacency(mtxA) ;
  /*
    ---------------------
    fill the Graph object
    ---------------------
  */
  Graph_init2(graph, 0, nvtx, 0, adjIVL->tsize, nvtx, adjIVL->tsize, 
              adjIVL, NULL, NULL) ;
  
  
  return(1) ;
}

/*--------------------------------------------------------------------*/



