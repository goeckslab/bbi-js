// Depends on jQuery, jQuery-ajax-native plugin, Underscore, jDataView.

// Adapted from JBrowse (https://github.com/GMOD/jbrowse)/src/JBrowse/Store/SeqFeature/BigWig.js
// which was adapted by Robert Buels from bigwig.js in the Dalliance Genome Explorer which is 
// copyright Thomas Down 2006-2010

/**
 * Provides data from a BigWig file.
 */
var BigWigDataProvider = function(options) {
    this.url = options.url;
};

// Add attributes/methods. 
_.extend(BigWigDataProvider.prototype, {
    BIG_WIG_MAGIC: -2003829722,
    BIG_BED_MAGIC: -2021002517,

    BIG_WIG_TYPE_GRAPH: 1,
    BIG_WIG_TYPE_VSTEP: 2,
    BIG_WIG_TYPE_FSTEP: 3,

    _littleEndian: true,

    chunkSizeLimit: 30000000, // 30mb

    /**
     * Load bigwig file.
     */
    _load: function() {
        var deferred = $.Deferred(),
            self = this;
        $.when(this._read(0, 512)).then(function(dataArraybuffer) {
            self._readHeader(dataArraybuffer);
            deferred.resolve();
        });
        return deferred;
    },

    /**
     * Read header.
     */
    _readHeader: function(dataArraybuffer) {
        var data = new jDataView(dataArraybuffer);

        // check magic numbers
        var magic = data.getInt32();
        if( magic !== this.BIG_WIG_MAGIC && magic !== this.BIG_BED_MAGIC ) {
            // try the other endianness if no magic
            this._littleEndian = true;
            data = this.newDataView( dataArraybuffer );
            magic = data.getInt32();
            if( magic !== this.BIG_WIG_MAGIC && magic !== this.BIG_BED_MAGIC) {
                console.error('Not a BigWig or BigBed file');
                deferred.reject('Not a BigWig or BigBed file');
                return;
            }
        }
        this.type = magic == this.BIG_BED_MAGIC ? 'bigbed' : 'bigwig';

        /*
        this.fileSize = bytes.fileSize;
        if( ! this.fileSize )
            console.warn("cannot get size of BigWig/BigBed file, widest zoom level not available");
        */

        this.version = data.getUint16();
        this.numZoomLevels = data.getUint16();
        this.chromTreeOffset = data.getUint64();
        this.unzoomedDataOffset = data.getUint64();
        this.unzoomedIndexOffset = data.getUint64();
        this.fieldCount = data.getUint16();
        this.definedFieldCount = data.getUint16();
        this.asOffset = data.getUint64();
        this.totalSummaryOffset = data.getUint64();
        this.uncompressBufSize = data.getUint32();

        // dlog('bigType: ' + this.type);
        // dlog('chromTree at: ' + this.chromTreeOffset);
        // dlog('uncompress: ' + this.uncompressBufSize);
        // dlog('data at: ' + this.unzoomedDataOffset);
        // dlog('index at: ' + this.unzoomedIndexOffset);
        // dlog('field count: ' + this.fieldCount);
        // dlog('defined count: ' + this.definedFieldCount);

        this.zoomLevels = [];
        for (var zl = 0; zl < this.numZoomLevels; ++zl) {
            var zlReduction = data.getUint32( 4*(zl*6 + 16) );
            var zlData = data.getUint64( 4*(zl*6 + 18) );
            var zlIndex = data.getUint64( 4*(zl*6 + 20) );

            //          dlog('zoom(' + zl + '): reduction=' + zlReduction + '; data=' + zlData + '; index=' + zlIndex);
            this.zoomLevels.push({reductionLevel: zlReduction, dataOffset: zlData, indexOffset: zlIndex});
        }

        // parse the totalSummary if present (summary of all data in the file)
        if( this.totalSummaryOffset ) {
            (function() {
                 var d = this.newDataView( dataArraybuffer, this.totalSummaryOffset.valueOf(), 512 - 304);
                 var s = {
                     basesCovered: d.getUint64(),
                     scoreMin: d.getFloat64(),
                     scoreMax: d.getFloat64(),
                     scoreSum: d.getFloat64(),
                     scoreSumSquares: d.getFloat64()
                 };
                 this._globalStats = s;
                 // rest of stats will be calculated on demand in getGlobalStats
             }).call(this);
        } else {
                console.warn("BigWig "+this.data.url+ " has no total summary data.");
        }

        // NOTE: this does not yet work.
        this._readChromTree(
            function() {
                //this._deferred.features.resolve({success: true});
                //this._deferred.stats.resolve({success: true});
            }
        );
    },

    _getGlobalStats: function( successCallback, errorCallback ) {
        var s = this._globalStats || {};

        // calc mean and standard deviation if necessary
        if( !( 'scoreMean' in s ))
            s.scoreMean = s.basesCovered ? s.scoreSum / s.basesCovered : 0;
        if( !( 'scoreStdDev' in s ))
            s.scoreStdDev = this._calcStdFromSums( s.scoreSum, s.scoreSumSquares, s.basesCovered );

        successCallback( s );
    },

     /**
      * Read from the bbi file, respecting the configured chunkSizeLimit.
      */
    _read: function(start, size, callback, errorcallback) {
        if(size > this.chunkSizeLimit) {
            // TODO: raise error.
        }
        else {
            this.read();
        }
    },

    newDataView: function( bytes, offset, length ) {
        return new jDataView( bytes, offset, length, this._littleEndian );
    },

    /**
     * @private
     */
    _readChromTree: function(callback, errorCallback) {
        var thisB = this;
        this.refsByNumber = {};
        this.refsByName = {};

        var udo = this.unzoomedDataOffset;
        while ((udo % 4) !== 0) {
            ++udo;
        }
        this._read( this.chromTreeOffset, udo - this.chromTreeOffset, function(bpt) {
            /*
           if( ! has('typed-arrays') ) {
               thisB._failAllDeferred( 'Web browser does not support typed arrays' );
               return;
           }
           */
           var data = thisB.newDataView( bpt );

           if( data.getUint32() !== 2026540177 )
               throw "parse error: not a Kent bPlusTree";
           var blockSize = data.getUint32();
           var keySize = data.getUint32();
           var valSize = data.getUint32();
           var itemCount = data.getUint64();
           var rootNodeOffset = 32;

           //dlog('blockSize=' + blockSize + '    keySize=' + keySize + '   valSize=' + valSize + '    itemCount=' + itemCount);

           var bptReadNode = function(offset) {
               if( offset >= bpt.length )
                   throw "reading beyond end of buffer";
               var isLeafNode = data.getUint8( offset );
               var cnt = data.getUint16( offset+2 );
               //dlog('ReadNode: ' + offset + '     type=' + isLeafNode + '   count=' + cnt);
               offset += 4;
               for (var n = 0; n < cnt; ++n) {
                   if( isLeafNode ) {
                       // parse leaf node
                       var key = '';
                       for (var ki = 0; ki < keySize; ++ki) {
                           var charCode = data.getUint8( offset++ );
                           if (charCode !== 0) {
                               key += String.fromCharCode(charCode);
                           }
                       }
                       var refId = data.getUint32( offset );
                       var refSize = data.getUint32( offset+4 );
                       offset += 8;

                       var refRec = { name: key, id: refId, length: refSize };

                       thisB.refsByName[ thisB.regularizeReferenceName(key) ] = refRec;
                       thisB.refsByNumber[refId] = refRec;
                   } else {
                       // parse index node
                       offset += keySize;
                       var childOffset = data.getUint64( offset );
                       offset += 8;
                       childOffset -= thisB.chromTreeOffset;
                       bptReadNode(childOffset);
                   }
               }
           };

           bptReadNode(rootNodeOffset);

           callback.call( thisB, thisB );
            }, errorCallback );
    },

    /**
     * Implement the JBrowse Browser.regularizeReferenceName
     * @param  {[type]} refname [description]
     * @return {[type]}         [description]
     */
    regularizeReferenceName: function(refname){
        /* TODO: Check if it is for save processing purpose or not, and if it is mendatory
        if( this.config.exactReferenceSequenceNames ){
            return refname;
        }
        */
        refname = refname.toLowerCase()
                         .replace(/^chro?m?(osome)?/,'chr')
                         .replace(/^co?n?ti?g/,'ctg')
                         .replace(/^scaff?o?l?d?/,'scaffold')
                         .replace(/^([a-z]*)0+/,'$1')
                         .replace(/^(\d+)$/, 'chr$1' );
        return refname;
    },

    /**
     * Interrogate whether a store has data for a given reference
     * sequence.  Calls the given callback with either true or false.
     *
     * Implemented as a binary interrogation because some stores are
     * smart enough to regularize reference sequence names, while
     * others are not.
     */
    hasRefSeq: function( seqName, callback, errorCallback ) {
        var thisB = this;
        seqName = thisB.browser.regularizeReferenceName( seqName );
        this._deferred.features.then(function() {
            callback( seqName in thisB.refsByName );
        }, errorCallback );
    },

    _getFeatures: function( query, featureCallback, endCallback, errorCallback ) {

        var chrName = this.regularizeReferenceName( query.ref );
        var min = query.start;
        var max = query.end;

        var v = query.basesPerSpan ? this.getView( 1/query.basesPerSpan ) :
                       query.scale ? this.getView( query.scale )          :
                                     this.getView( 1 );

        if( !v ) {
            endCallback();
            return;
        }

        v.readWigData( chrName, min, max, dojo.hitch( this, function( features ) {
            array.forEach( features || [], featureCallback );
            endCallback();
        }), errorCallback );
    },

    readWigData: function(chrName, min, max, callback, errorCallback ) {
        // console.log( 'reading wig data from '+chrName+':'+min+'..'+max);
        var chr = this.bwg.refsByName[chrName];
        if ( ! chr ) {
            // Not an error because some .bwgs won't have data for all chromosomes.

            // dlog("Couldn't find chr " + chrName);
            // dlog('Chroms=' + miniJSONify(this.bwg.refsByName));
            callback([]);
        } else {
            this.readWigDataById( chr.id, min, max, callback, errorCallback );
        }
    },

    readWigDataById: function(chr, min, max, callback, errorCallback ) {
        // TODO: this was representing Window
        if( !this.cirHeader ) {
            var readCallback = lang.hitch( this, 'readWigDataById', chr, min, max, callback, errorCallback );
            if( this.cirHeaderLoading ) {
                this.cirHeaderLoading.push( readCallback );
            }
            else {
                this.cirHeaderLoading = [ readCallback ];
                // dlog('No CIR yet, fetching');
                this._read( this.cirTreeOffset, 48, lang.hitch( this, function(result) {
                                this.cirHeader = result;
                                this.cirBlockSize = this.newDataView( result, 4, 4 ).getUint32();
                                array.forEach( this.cirHeaderLoading, function(c) { c(); });
                                delete this.cirHeaderLoading;
                            }), errorCallback );
            }
            return;
        }

        //dlog('_readWigDataById', chr, min, max, callback);

        var worker = new RequestWorker( this, chr, min, max, callback, errorCallback );
        worker.cirFobRecur([this.cirTreeOffset + 48], 1);
    },

    getUnzoomedView: function() {
        if (!this.unzoomedView) {
            var cirLen = 4000;
            var nzl = this.zoomLevels[0];
            if (nzl) {
                cirLen = this.zoomLevels[0].dataOffset - this.unzoomedIndexOffset;
            }
            this.unzoomedView = new Window( this, this.unzoomedIndexOffset, cirLen, false );
        }
        return this.unzoomedView;
    },

    getView: function( scale ) {
        if( ! this.zoomLevels || ! this.zoomLevels.length )
            return null;

        if( !this._viewCache || this._viewCache.scale != scale ) {
            this._viewCache = {
                scale: scale,
                view: this._getView( scale )
            };
        }
        return this._viewCache.view;
    },

    _getView: function( scale ) {
        var basesPerPx = 1/scale;
        //console.log('getting view for '+basesPerSpan+' bases per span');
        var maxLevel = this.zoomLevels.length;
        if( ! this.fileSize ) // if we don't know the file size, we can't fetch the highest zoom level :-(
            maxLevel--;
        for( var i = maxLevel; i > 0; i-- ) {
            var zh = this.zoomLevels[i];
            if( zh && zh.reductionLevel <= 2*basesPerPx ) {
                var indexLength = i < this.zoomLevels.length - 1 ?
                    this.zoomLevels[i + 1].dataOffset - zh.indexOffset
                    : this.fileSize - 4 - zh.indexOffset;
                //console.log( 'using zoom level '+i);
                return new Window( this, zh.indexOffset, indexLength, true );
            }
        }
        //console.log( 'using unzoomed level');
        return this.getUnzoomedView();
    }
});

