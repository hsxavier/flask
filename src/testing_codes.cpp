// Test allocation of gsl_matrix array.

  CovByl=GSLMatrixArray(3, 2, 2);
  // [i*L->size1+j]
  l=0;
  CovByl[l]->data[0*2+0]=10*(l+1)+1; CovByl[l]->data[0*2+1]=10*(l+1)+2;
  CovByl[l]->data[1*2+0]=10*(l+1)+3; CovByl[l]->data[1*2+1]=10*(l+1)+4;
  l=1;
  CovByl[l]->data[0*2+0]=10*(l+1)+1; CovByl[l]->data[0*2+1]=10*(l+1)+2;
  CovByl[l]->data[1*2+0]=10*(l+1)+3; CovByl[l]->data[1*2+1]=10*(l+1)+4;
  l=2;
  CovByl[l]->data[0*2+0]=10*(l+1)+1; CovByl[l]->data[0*2+1]=10*(l+1)+2;
  CovByl[l]->data[1*2+0]=10*(l+1)+3; CovByl[l]->data[1*2+1]=10*(l+1)+4;
  
  for (l=0; l<3; l++) {
    for (i=0; i<2; i++) {
      for (j=0; j<2; j++) printf("%g ",CovByl[l]->data[i*2+j]);
      printf("\n");
    }
    printf("\n");
  }
  return 0;
