from proxigenomics_toolkit.misc_utils import *
from proxigenomics_toolkit.io_utils import *
from proxigenomics_toolkit.contact_map import *
from proxigenomics_toolkit.exceptions import ApplicationException
import logging
import sys

__version__ = '0.1'

if __name__ == '__main__':
    import argparse

    def mk_version():
        return 'scaffold3C v{}'.format(__version__)

    def out_name(base, suffix):
        return '{}{}'.format(base, suffix)

    parser = argparse.ArgumentParser(description='Create a 3C fragment map from a BAM file')

    parser.add_argument('-V', '--version', default=False, action='store_true', help='Show the application version')
    parser.add_argument('-s', '--seed', default=None, type=int, help='Random integer seed')
    parser.add_argument('-v', '--verbose', default=False, action='store_true', help='Verbose output')
    parser.add_argument('--clobber', default=False, action='store_true', help='Clobber existing files')
    parser.add_argument('--log', help='Log file path [OUTDIR/scaffold3C.log]')
    parser.add_argument('--max-image', type=int, default=4000, help='Maximum image size for plots [4000]')

    parser.add_argument('--min-reflen', type=int, default=1000, help='Minimum acceptable reference length [1000]')
    parser.add_argument('--min-signal', type=int, default=5, help='Minimum acceptable trans signal [5]')
    parser.add_argument('--min-size', type=int, default=5, help='Minimum cluster size for ordering [5]')
    parser.add_argument('--min-extent', type=int, default=50000,
                        help='Minimum cluster extent (kb) for ordering [50000]')
    parser.add_argument('--min-ordlen', default=2000,
                        help='Minimum length of sequence to use in ordering [2000]')

    parser.add_argument('--dist-method', choices=['inverse', 'neglog'], default='inverse',
                        help='Distance method for ordering [inverse]')
    parser.add_argument('--only-large', default=False, action='store_true',
                        help='Only write FASTA for clusters longer than min_extent')
    parser.add_argument('--skip-plotting', default=False, action='store_true',
                        help='Skip plotting the contact map')
    parser.add_argument('--fasta', help='Alternative location of source FASTA from that supplied during clustering')
    parser.add_argument('MAP', help='Contact map')
    parser.add_argument('CLUSTERING', help='Clustering solution')
    parser.add_argument('OUTDIR', help='Output directory')

    args = parser.parse_args()

    if args.version:
        print mk_version()
        sys.exit(0)

    try:
        make_dir(args.OUTDIR, args.clobber)
    except IOError as e:
        print 'Error: {}'.format(e.message)
        sys.exit(1)

    logging.captureWarnings(True)
    logger = logging.getLogger('main')

    # root log listens to everything
    root = logging.getLogger('')
    root.setLevel(logging.DEBUG)

    # log message format
    formatter = logging.Formatter(fmt='%(levelname)-8s | %(asctime)s | %(name)7s | %(message)s')

    # Runtime console listens to INFO by default
    ch = logging.StreamHandler()
    if args.verbose:
        ch.setLevel(logging.DEBUG)
    else:
        ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    root.addHandler(ch)

    # File log listens to all levels from root
    if args.log is not None:
        log_path = args.log
    else:
        log_path = os.path.join(args.OUTDIR, 'scaffold3C.log')
    fh = logging.FileHandler(log_path, mode='a')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    root.addHandler(fh)

    # Add some environmental details
    logger.debug(mk_version())
    logger.debug(sys.version.replace('\n', ' '))
    logger.debug('Command line: {}'.format(' '.join(sys.argv)))

    try:

        make_dir(args.OUTDIR, exist_ok=True)

        if not args.seed:
            args.seed = make_random_seed()
            logger.info('Generated random seed: {}'.format(args.seed))
        else:
            logger.info("User set random seed: {}".format(args.seed))

        # Load a pre-existing serialized contact map and clustering solution
        logger.info('Loading existing contact map from: {}'.format(args.MAP))
        cm = load_object(args.MAP)

        logger.info('Loading existing clustering solution from: {}'.format(args.CLUSTERING))
        clustering = load_object(args.CLUSTERING)

        cm.min_extent = args.min_extent
        # update the mask if the user has specified a threshold, which might be the same as before
        if args.min_signal != cm.min_sig or args.min_reflen != cm.min_len:
            # pedantically set these and pass to method just in-case of logic oversight
            cm.min_len = args.min_reflen
            cm.min_sig = args.min_signal
            cm.set_primary_acceptance_mask(min_sig=args.min_signal, min_len=args.min_reflen, update=True)

        if cm.is_empty():
            logger.info('Stopping as the map is empty')
            sys.exit(1)

        # order
        order_clusters(cm, clustering, seed=args.seed, min_len=args.min_ordlen, min_size=args.min_size,
                       dist_method=args.dist_method, work_dir=args.OUTDIR)

        # serialize full clustering object again
        logger.info('Saving ordered clustering instance')
        save_object(os.path.join(args.OUTDIR, 'clustering_ordered.p'), clustering)

        # write per-cluster fasta files, also separate ordered fasta if ordering performed
        write_fasta(cm, args.OUTDIR, clustering, source_fasta=args.fasta, clobber=True, only_large=args.only_large)

        if not args.skip_plotting:

            # just the clusters and contigs which were ordered
            logger.info('Plotting ordered clusters')
            plot_clusters(cm, os.path.join(args.OUTDIR, 'cluster_scaffolded_plot.png'), clustering,
                          max_image_size=args.max_image, ordered_only=True, simple=False, permute=True)

    except ApplicationException as ex:
        import sys
        logger.error(ex.message)
        sys.exit(1)
