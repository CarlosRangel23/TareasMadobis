### GRUPO 4

## English -- Spanish at the bottom
This script automatically generates, for a limited number of parameters and conditions, the .narrowPeak and .bed files for a chipseq sample. It is important to note that the file does not work if we do not provide any input file, just as it would not work if, in case there is more than one input, the number of inputs differs from the number of samples treated. In addition, .bam and .bw files are also generated for each sample, including the control, which are not deleted in case there is interest in analyzing them in software such as the IGV.

Before running the code, the following experimental design must be available:

* annotation
    * reference_genome.gtf
* genome
    * genome.fa
* results
* samples
    * chip01
    * chip02
    * ...
    * input01
    * input02
    * ...
* scripts
    * onescript_forall.sh
    * sampledir
    * otherdir
    * sampleSRA

How should the sampledir, otherdir and sampleSRA files be organized?
- sampledir: must contain the absolute paths to the folder where you want the samples to be downloaded and the .bam .bam.bai and .bw files to be generated. It is important that in the sampledir file one absolute path is written per line, as shown in the example. On the other hand, it is also important to order the routes correctly, that is, if we have several inputs, maintain the same order between the treated samples and the inputs. Thus, when peak calling is made, the macs2 program can be used sequentially, using its corresponding input for each sample. Finally, it is important to place the route of the treated samples first and then the route of the inputs.
- otherdir: it is another file where we will store the file paths. The first absolute path to be written will correspond to genome, the second to the folder where sampledir, otherdir and sampleSRA are located. Finally, the third route will be the one that indicates where we want the results to be saved (results in our case).
- sampleSRA: the accession number of each sample, including the inputs. Each line corresponds to a sample and it is important to follow the same order as the one we have followed in sampledir.

An example is shown in the [test](https://github.com/CarlosRangel23/TareasMadobis/tree/main/OmicsSubject/chipseq/test) folder.

## Espa??ol (Spanish)
Este script genera autom??ticamente, para un n??mero limitado de par??metros y condiciones, los archivos .narrowPeak y .bed para una muestra de chipseq. Es importante tener en cuenta que el archivo no funciona si no proporcionamos ning??n archivo input, al igual que tampoco funcionar??a si, en caso de que haya m??s de un input, el n??mero de inputs difiera del n??mero de muestras tratadas. Adem??s, tambi??n se genera para cada muestra, incluida el control, archivos .bam y .bw, que no son eliminados por si existe inter??s en analizarlos en softwares como el IGV.

Antes de correr el c??digo, se debe disponer del siguiente dise??o experimental:

* annotation
    * reference_genome.gtf
* genome
    * genome.fa
* results
* samples
    * chip01
    * chip02
    * ...
    * input01
    * input02
    * ...
* scripts
    * onescript_forall.sh
    * sampledir
    * otherdir
    * sampleSRA
     
??C??mo deben organizarse los archivos sampledir, otherdir y sampleSRA?
- sampledir: debe contener las rutas absolutas a la carpeta donde quieres que se descarguen las muestras y se generen los archivos .bam .bam.bai y .bw. Es importante que en el fichero sampledir se escriba una ruta absoluta por l??nea, como se muestra en el ejemplo. Por otra parte, es importante tambi??n ordenar las rutas correctamente, es decir, si disponemos de varios inputs, mantener el mismo orden entre las muestras tratadas y los inputs. As??, cuando se haga el peak calling se podr?? utilizar secuencialmente el programa macs2 utilizando para cada muestra su correspondiente input. Finalmente, es importante colocar primero la ruta de las muestras tratdas y posteriormente la ruta de los inputs.
- otherdir: es otro fichero donde guardaremos las rutas de archivos. La primera ruta absoluta que se escribir?? corresponder?? a genome, la segunda al folder donde se encuentren sampledir, otherdir y sampleSRA. Por ??ltimo, la tercera ruta ser?? aquella que indique donde queremos que se guarden los resultados (results en nuestro caso).
- sampleSRA: los accession number de cada muestra, incluyendo los inputs. Cada l??nea corresponde con una muestra y es importante que se siga el mismo orden que aquel que hemos seguido en sampledir.

Un ejemplo se muestra en la carpeta [test](https://github.com/CarlosRangel23/TareasMadobis/tree/main/OmicsSubject/chipseq/test).
