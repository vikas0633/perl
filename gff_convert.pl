### perl script to interchange between various gff format

use Bio::FeatureIO;

#read from a file
$in  = Bio::FeatureIO->new(-file => "lotus_pasa.validated_transcripts.gff3" , -format => 'GFF', -version => 3);

#read from a filehandle
$in  = Bio::FeatureIO->new(-fh => \*GFF , -format => 'GFF');

#read features already attached to a sequence
my $feat = Bio::FeatureIO->new(-seq => $seq , -format => 'features');

#read new features for existing sequence
my $seq = Bio::FeatureIO->new(-seq => $seq , -format => 'Das');

#write out features
$out = Bio::FeatureIO->new(-file    => ">lotus_pasa.validated_transcripts.gff" ,
						 -format  => 'GFF' , );

while ( my $feature = $in->next_feature() ) {
$out->write_feature($feature);
}