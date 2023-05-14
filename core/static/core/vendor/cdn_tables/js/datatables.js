/*
 * This combined file was created by the DataTables downloader builder:
 *   https://datatables.net/download
 *
 * To rebuild or modify this file with the latest versions of the included
 * software please visit:
 *   https://datatables.net/download/#bs5/dt-1.13.4/b-2.3.6/sb-1.4.2/sp-2.1.2
 *
 * Included libraries:
 *   DataTables 1.13.4, Buttons 2.3.6, SearchBuilder 1.4.2, SearchPanes 2.1.2
 */

/*! DataTables 1.13.4
 * ©2008-2023 SpryMedia Ltd - datatables.net/license
 */

/**
 * @summary     DataTables
 * @description Paginate, search and order HTML tables
 * @version     1.13.4
 * @author      SpryMedia Ltd
 * @contact     www.datatables.net
 * @copyright   SpryMedia Ltd.
 *
 * This source file is free software, available under the following license:
 *   MIT license - http://datatables.net/license
 *
 * This source file is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the license files for details.
 *
 * For details please refer to: http://www.datatables.net
 */

/*jslint evil: true, undef: true, browser: true */
/*globals $,require,jQuery,define,_selector_run,_selector_opts,_selector_first,_selector_row_indexes,_ext,_Api,_api_register,_api_registerPlural,_re_new_lines,_re_html,_re_formatted_numeric,_re_escape_regex,_empty,_intVal,_numToDecimal,_isNumber,_isHtml,_htmlNumeric,_pluck,_pluck_order,_range,_stripHtml,_unique,_fnBuildAjax,_fnAjaxUpdate,_fnAjaxParameters,_fnAjaxUpdateDraw,_fnAjaxDataSrc,_fnAddColumn,_fnColumnOptions,_fnAdjustColumnSizing,_fnVisibleToColumnIndex,_fnColumnIndexToVisible,_fnVisbleColumns,_fnGetColumns,_fnColumnTypes,_fnApplyColumnDefs,_fnHungarianMap,_fnCamelToHungarian,_fnLanguageCompat,_fnBrowserDetect,_fnAddData,_fnAddTr,_fnNodeToDataIndex,_fnNodeToColumnIndex,_fnGetCellData,_fnSetCellData,_fnSplitObjNotation,_fnGetObjectDataFn,_fnSetObjectDataFn,_fnGetDataMaster,_fnClearTable,_fnDeleteIndex,_fnInvalidate,_fnGetRowElements,_fnCreateTr,_fnBuildHead,_fnDrawHead,_fnDraw,_fnReDraw,_fnAddOptionsHtml,_fnDetectHeader,_fnGetUniqueThs,_fnFeatureHtmlFilter,_fnFilterComplete,_fnFilterCustom,_fnFilterColumn,_fnFilter,_fnFilterCreateSearch,_fnEscapeRegex,_fnFilterData,_fnFeatureHtmlInfo,_fnUpdateInfo,_fnInfoMacros,_fnInitialise,_fnInitComplete,_fnLengthChange,_fnFeatureHtmlLength,_fnFeatureHtmlPaginate,_fnPageChange,_fnFeatureHtmlProcessing,_fnProcessingDisplay,_fnFeatureHtmlTable,_fnScrollDraw,_fnApplyToChildren,_fnCalculateColumnWidths,_fnThrottle,_fnConvertToWidth,_fnGetWidestNode,_fnGetMaxLenString,_fnStringToCss,_fnSortFlatten,_fnSort,_fnSortAria,_fnSortListener,_fnSortAttachListener,_fnSortingClasses,_fnSortData,_fnSaveState,_fnLoadState,_fnSettingsFromNode,_fnLog,_fnMap,_fnBindAction,_fnCallbackReg,_fnCallbackFire,_fnLengthOverflow,_fnRenderer,_fnDataSource,_fnRowAttributes*/

(function( factory ) {
	"use strict";

	if ( typeof define === 'function' && define.amd ) {
		// AMD
		define( ['jquery'], function ( $ ) {
			return factory( $, window, document );
		} );
	}
	else if ( typeof exports === 'object' ) {
		// CommonJS
		// jQuery's factory checks for a global window - if it isn't present then it
		// returns a factory function that expects the window object
		var jq = require('jquery');

		if (typeof window !== 'undefined') {
			module.exports = function (root, $) {
				if ( ! root ) {
					// CommonJS environments without a window global must pass a
					// root. This will give an error otherwise
					root = window;
				}

				if ( ! $ ) {
					$ = jq( root );
				}

				return factory( $, root, root.document );
			};
		}
		else {
			return factory( jq, window, window.document );
		}
	}
	else {
		// Browser
		window.DataTable = factory( jQuery, window, document );
	}
}
(function( $, window, document, undefined ) {
	"use strict";

	
	var DataTable = function ( selector, options )
	{
		// Check if called with a window or jQuery object for DOM less applications
		// This is for backwards compatibility
		if (DataTable.factory(selector, options)) {
			return DataTable;
		}
	
		// When creating with `new`, create a new DataTable, returning the API instance
		if (this instanceof DataTable) {
			return $(selector).DataTable(options);
		}
		else {
			// Argument switching
			options = selector;
		}
	
		/**
		 * Perform a jQuery selector action on the table's TR elements (from the tbody) and
		 * return the resulting jQuery object.
		 *  @param {string|node|jQuery} sSelector jQuery selector or node collection to act on
		 *  @param {object} [oOpts] Optional parameters for modifying the rows to be included
		 *  @param {string} [oOpts.filter=none] Select TR elements that meet the current filter
		 *    criterion ("applied") or all TR elements (i.e. no filter).
		 *  @param {string} [oOpts.order=current] Order of the TR elements in the processed array.
		 *    Can be either 'current', whereby the current sorting of the table is used, or
		 *    'original' whereby the original order the data was read into the table is used.
		 *  @param {string} [oOpts.page=all] Limit the selection to the currently displayed page
		 *    ("current") or not ("all"). If 'current' is given, then order is assumed to be
		 *    'current' and filter is 'applied', regardless of what they might be given as.
		 *  @returns {object} jQuery object, filtered by the given selector.
		 *  @dtopt API
		 *  @deprecated Since v1.10
		 *
		 *  @example
		 *    $(document).ready(function() {
		 *      var oTable = $('#example').dataTable();
		 *
		 *      // Highlight every second row
		 *      oTable.$('tr:odd').css('backgroundColor', 'blue');
		 *    } );
		 *
		 *  @example
		 *    $(document).ready(function() {
		 *      var oTable = $('#example').dataTable();
		 *
		 *      // Filter to rows with 'Webkit' in them, add a background colour and then
		 *      // remove the filter, thus highlighting the 'Webkit' rows only.
		 *      oTable.fnFilter('Webkit');
		 *      oTable.$('tr', {"search": "applied"}).css('backgroundColor', 'blue');
		 *      oTable.fnFilter('');
		 *    } );
		 */
		this.$ = function ( sSelector, oOpts )
		{
			return this.api(true).$( sSelector, oOpts );
		};
		
		
		/**
		 * Almost identical to $ in operation, but in this case returns the data for the matched
		 * rows - as such, the jQuery selector used should match TR row nodes or TD/TH cell nodes
		 * rather than any descendants, so the data can be obtained for the row/cell. If matching
		 * rows are found, the data returned is the original data array/object that was used to
		 * create the row (or a generated array if from a DOM source).
		 *
		 * This method is often useful in-combination with $ where both functions are given the
		 * same parameters and the array indexes will match identically.
		 *  @param {string|node|jQuery} sSelector jQuery selector or node collection to act on
		 *  @param {object} [oOpts] Optional parameters for modifying the rows to be included
		 *  @param {string} [oOpts.filter=none] Select elements that meet the current filter
		 *    criterion ("applied") or all elements (i.e. no filter).
		 *  @param {string} [oOpts.order=current] Order of the data in the processed array.
		 *    Can be either 'current', whereby the current sorting of the table is used, or
		 *    'original' whereby the original order the data was read into the table is used.
		 *  @param {string} [oOpts.page=all] Limit the selection to the currently displayed page
		 *    ("current") or not ("all"). If 'current' is given, then order is assumed to be
		 *    'current' and filter is 'applied', regardless of what they might be given as.
		 *  @returns {array} Data for the matched elements. If any elements, as a result of the
		 *    selector, were not TR, TD or TH elements in the DataTable, they will have a null
		 *    entry in the array.
		 *  @dtopt API
		 *  @deprecated Since v1.10
		 *
		 *  @example
		 *    $(document).ready(function() {
		 *      var oTable = $('#example').dataTable();
		 *
		 *      // Get the data from the first row in the table
		 *      var data = oTable._('tr:first');
		 *
		 *      // Do something useful with the data
		 *      alert( "First cell is: "+data[0] );
		 *    } );
		 *
		 *  @example
		 *    $(document).ready(function() {
		 *      var oTable = $('#example').dataTable();
		 *
		 *      // Filter to 'Webkit' and get all data for
		 *      oTable.fnFilter('Webkit');
		 *      var data = oTable._('tr', {"search": "applied"});
		 *
		 *      // Do something with the data
		 *      alert( data.length+" rows matched the search" );
		 *    } );
		 */
		this._ = function ( sSelector, oOpts )
		{
			return this.api(true).rows( sSelector, oOpts ).data();
		};
		
		
		/**
		 * Create a DataTables Api instance, with the currently selected tables for
		 * the Api's context.
		 * @param {boolean} [traditional=false] Set the API instance's context to be
		 *   only the table referred to by the `DataTable.ext.iApiIndex` option, as was
		 *   used in the API presented by DataTables 1.9- (i.e. the traditional mode),
		 *   or if all tables captured in the jQuery object should be used.
		 * @return {DataTables.Api}
		 */
		this.api = function ( traditional )
		{
			return traditional ?
				new _Api(
					_fnSettingsFromNode( this[ _ext.iApiIndex ] )
				) :
				new _Api( this );
		};
		
		
		/**
		 * Add a single new row or multiple rows of data to the table. Please note
		 * that this is suitable for client-side processing only - if you are using
		 * server-side processing (i.e. "bServerSide": true), then to add data, you
		 * must add it to the data source, i.e. the server-side, through an Ajax call.
		 *  @param {array|object} data The data to be added to the table. This can be:
		 *    <ul>
		 *      <li>1D array of data - add a single row with the data provided</li>
		 *      <li>2D array of arrays - add multiple rows in a single call</li>
		 *      <li>object - data object when using <i>mData</i></li>
		 *      <li>array of objects - multiple data objects when using <i>mData</i></li>
		 *    </ul>
		 *  @param {bool} [redraw=true] redraw the table or not
		 *  @returns {array} An array of integers, representing the list of indexes in
		 *    <i>aoData</i> ({@link DataTable.models.oSettings}) that have been added to
		 *    the table.
		 *  @dtopt API
		 *  @deprecated Since v1.10
		 *
		 *  @example
		 *    // Global var for counter
		 *    var giCount = 2;
		 *
		 *    $(document).ready(function() {
		 *      $('#example').dataTable();
		 *    } );
		 *
		 *    function fnClickAddRow() {
		 *      $('#example').dataTable().fnAddData( [
		 *        giCount+".1",
		 *        giCount+".2",
		 *        giCount+".3",
		 *        giCount+".4" ]
		 *      );
		 *
		 *      giCount++;
		 *    }
		 */
		this.fnAddData = function( data, redraw )
		{
			var api = this.api( true );
		
			/* Check if we want to add multiple rows or not */
			var rows = Array.isArray(data) && ( Array.isArray(data[0]) || $.isPlainObject(data[0]) ) ?
				api.rows.add( data ) :
				api.row.add( data );
		
			if ( redraw === undefined || redraw ) {
				api.draw();
			}
		
			return rows.flatten().toArray();
		};
		
		
		/**
		 * This function will make DataTables recalculate the column sizes, based on the data
		 * contained in the table and the sizes applied to the columns (in the DOM, CSS or
		 * through the sWidth parameter). This can be useful when the width of the table's
		 * parent element changes (for example a window resize).
		 *  @param {boolean} [bRedraw=true] Redraw the table or not, you will typically want to
		 *  @dtopt API
		 *  @deprecated Since v1.10
		 *
		 *  @example
		 *    $(document).ready(function() {
		 *      var oTable = $('#example').dataTable( {
		 *        "sScrollY": "200px",
		 *        "bPaginate": false
		 *      } );
		 *
		 *      $(window).on('resize', function () {
		 *        oTable.fnAdjustColumnSizing();
		 *      } );
		 *    } );
		 */
		this.fnAdjustColumnSizing = function ( bRedraw )
		{
			var api = this.api( true ).columns.adjust();
			var settings = api.settings()[0];
			var scroll = settings.oScroll;
		
			if ( bRedraw === undefined || bRedraw ) {
				api.draw( false );
			}
			else if ( scroll.sX !== "" || scroll.sY !== "" ) {
				/* If not redrawing, but scrolling, we want to apply the new column sizes anyway */
				_fnScrollDraw( settings );
			}
		};
		
		
		/**
		 * Quickly and simply clear a table
		 *  @param {bool} [bRedraw=true] redraw the table or not
		 *  @dtopt API
		 *  @deprecated Since v1.10
		 *
		 *  @example
		 *    $(document).ready(function() {
		 *      var oTable = $('#example').dataTable();
		 *
		 *      // Immediately 'nuke' the current rows (perhaps waiting for an Ajax callback...)
		 *      oTable.fnClearTable();
		 *    } );
		 */
		this.fnClearTable = function( bRedraw )
		{
			var api = this.api( true ).clear();
		
			if ( bRedraw === undefined || bRedraw ) {
				api.draw();
			}
		};
		
		
		/**
		 * The exact opposite of 'opening' a row, this function will close any rows which
		 * are currently 'open'.
		 *  @param {node} nTr the table row to 'close'
		 *  @returns {int} 0 on success, or 1 if failed (can't find the row)
		 *  @dtopt API
		 *  @deprecated Since v1.10
		 *
		 *  @example
		 *    $(document).ready(function() {
		 *      var oTable;
		 *
		 *      // 'open' an information row when a row is clicked on
		 *      $('#example tbody tr').click( function () {
		 *        if ( oTable.fnIsOpen(this) ) {
		 *          oTable.fnClose( this );
		 *        } else {
		 *          oTable.fnOpen( this, "Temporary row opened", "info_row" );
		 *        }
		 *      } );
		 *
		 *      oTable = $('#example').dataTable();
		 *    } );
		 */
		this.fnClose = function( nTr )
		{
			this.api( true ).row( nTr ).child.hide();
		};
		
		
		/**
		 * Remove a row for the table
		 *  @param {mixed} target The index of the row from aoData to be deleted, or
		 *    the TR element you want to delete
		 *  @param {function|null} [callBack] Callback function
		 *  @param {bool} [redraw=true] Redraw the table or not
		 *  @returns {array} The row that was deleted
		 *  @dtopt API
		 *  @deprecated Since v1.10
		 *
		 *  @example
		 *    $(document).ready(function() {
		 *      var oTable = $('#example').dataTable();
		 *
		 *      // Immediately remove the first row
		 *      oTable.fnDeleteRow( 0 );
		 *    } );
		 */
		this.fnDeleteRow = function( target, callback, redraw )
		{
			var api = this.api( true );
			var rows = api.rows( target );
			var settings = rows.settings()[0];
			var data = settings.aoData[ rows[0][0] ];
		
			rows.remove();
		
			if ( callback ) {
				callback.call( this, settings, data );
			}
		
			if ( redraw === undefined || redraw ) {
				api.draw();
			}
		
			return data;
		};
		
		
		/**
		 * Restore the table to it's original state in the DOM by removing all of DataTables
		 * enhancements, alterations to the DOM structure of the table and event listeners.
		 *  @param {boolean} [remove=false] Completely remove the table from the DOM
		 *  @dtopt API
		 *  @deprecated Since v1.10
		 *
		 *  @example
		 *    $(document).ready(function() {
		 *      // This example is fairly pointless in reality, but shows how fnDestroy can be used
		 *      var oTable = $('#example').dataTable();
		 *      oTable.fnDestroy();
		 *    } );
		 */
		this.fnDestroy = function ( remove )
		{
			this.api( true ).destroy( remove );
		};
		
		
		/**
		 * Redraw the table
		 *  @param {bool} [complete=true] Re-filter and resort (if enabled) the table before the draw.
		 *  @dtopt API
		 *  @deprecated Since v1.10
		 *
		 *  @example
		 *    $(document).ready(function() {
		 *      var oTable = $('#example').dataTable();
		 *
		 *      // Re-draw the table - you wouldn't want to do it here, but it's an example :-)
		 *      oTable.fnDraw();
		 *    } );
		 */
		this.fnDraw = function( complete )
		{
			// Note that this isn't an exact match to the old call to _fnDraw - it takes
			// into account the new data, but can hold position.
			this.api( true ).draw( complete );
		};
		
		
		/**
		 * Filter the input based on data
		 *  @param {string} sInput String to filter the table on
		 *  @param {int|null} [iColumn] Column to limit filtering to
		 *  @param {bool} [bRegex=false] Treat as regular expression or not
		 *  @param {bool} [bSmart=true] Perform smart filtering or not
		 *  @param {bool} [bShowGlobal=true] Show the input global filter in it's input box(es)
		 *  @param {bool} [bCaseInsensitive=true] Do case-insensitive matching (true) or not (false)
		 *  @dtopt API
		 *  @deprecated Since v1.10
		 *
		 *  @example
		 *    $(document).ready(function() {
		 *      var oTable = $('#example').dataTable();
		 *
		 *      // Sometime later - filter...
		 *      oTable.fnFilter( 'test string' );
		 *    } );
		 */
		this.fnFilter = function( sInput, iColumn, bRegex, bSmart, bShowGlobal, bCaseInsensitive )
		{
			var api = this.api( true );
		
			if ( iColumn === null || iColumn === undefined ) {
				api.search( sInput, bRegex, bSmart, bCaseInsensitive );
			}
			else {
				api.column( iColumn ).search( sInput, bRegex, bSmart, bCaseInsensitive );
			}
		
			api.draw();
		};
		
		
		/**
		 * Get the data for the whole table, an individual row or an individual cell based on the
		 * provided parameters.
		 *  @param {int|node} [src] A TR row node, TD/TH cell node or an integer. If given as
		 *    a TR node then the data source for the whole row will be returned. If given as a
		 *    TD/TH cell node then iCol will be automatically calculated and the data for the
		 *    cell returned. If given as an integer, then this is treated as the aoData internal
		 *    data index for the row (see fnGetPosition) and the data for that row used.
		 *  @param {int} [col] Optional column index that you want the data of.
		 *  @returns {array|object|string} If mRow is undefined, then the data for all rows is
		 *    returned. If mRow is defined, just data for that row, and is iCol is
		 *    defined, only data for the designated cell is returned.
		 *  @dtopt API
		 *  @deprecated Since v1.10
		 *
		 *  @example
		 *    // Row data
		 *    $(document).ready(function() {
		 *      oTable = $('#example').dataTable();
		 *
		 *      oTable.$('tr').click( function () {
		 *        var data = oTable.fnGetData( this );
		 *        // ... do something with the array / object of data for the row
		 *      } );
		 *    } );
		 *
		 *  @example
		 *    // Individual cell data
		 *    $(document).ready(function() {
		 *      oTable = $('#example').dataTable();
		 *
		 *      oTable.$('td').click( function () {
		 *        var sData = oTable.fnGetData( this );
		 *        alert( 'The cell clicked on had the value of '+sData );
		 *      } );
		 *    } );
		 */
		this.fnGetData = function( src, col )
		{
			var api = this.api( true );
		
			if ( src !== undefined ) {
				var type = src.nodeName ? src.nodeName.toLowerCase() : '';
		
				return col !== undefined || type == 'td' || type == 'th' ?
					api.cell( src, col ).data() :
					api.row( src ).data() || null;
			}
		
			return api.data().toArray();
		};
		
		
		/**
		 * Get an array of the TR nodes that are used in the table's body. Note that you will
		 * typically want to use the '$' API method in preference to this as it is more
		 * flexible.
		 *  @param {int} [iRow] Optional row index for the TR element you want
		 *  @returns {array|node} If iRow is undefined, returns an array of all TR elements
		 *    in the table's body, or iRow is defined, just the TR element requested.
		 *  @dtopt API
		 *  @deprecated Since v1.10
		 *
		 *  @example
		 *    $(document).ready(function() {
		 *      var oTable = $('#example').dataTable();
		 *
		 *      // Get the nodes from the table
		 *      var nNodes = oTable.fnGetNodes( );
		 *    } );
		 */
		this.fnGetNodes = function( iRow )
		{
			var api = this.api( true );
		
			return iRow !== undefined ?
				api.row( iRow ).node() :
				api.rows().nodes().flatten().toArray();
		};
		
		
		/**
		 * Get the array indexes of a particular cell from it's DOM element
		 * and column index including hidden columns
		 *  @param {node} node this can either be a TR, TD or TH in the table's body
		 *  @returns {int} If nNode is given as a TR, then a single index is returned, or
		 *    if given as a cell, an array of [row index, column index (visible),
		 *    column index (all)] is given.
		 *  @dtopt API
		 *  @deprecated Since v1.10
		 *
		 *  @example
		 *    $(document).ready(function() {
		 *      $('#example tbody td').click( function () {
		 *        // Get the position of the current data from the node
		 *        var aPos = oTable.fnGetPosition( this );
		 *
		 *        // Get the data array for this row
		 *        var aData = oTable.fnGetData( aPos[0] );
		 *
		 *        // Update the data array and return the value
		 *        aData[ aPos[1] ] = 'clicked';
		 *        this.innerHTML = 'clicked';
		 *      } );
		 *
		 *      // Init DataTables
		 *      oTable = $('#example').dataTable();
		 *    } );
		 */
		this.fnGetPosition = function( node )
		{
			var api = this.api( true );
			var nodeName = node.nodeName.toUpperCase();
		
			if ( nodeName == 'TR' ) {
				return api.row( node ).index();
			}
			else if ( nodeName == 'TD' || nodeName == 'TH' ) {
				var cell = api.cell( node ).index();
		
				return [
					cell.row,
					cell.columnVisible,
					cell.column
				];
			}
			return null;
		};
		
		
		/**
		 * Check to see if a row is 'open' or not.
		 *  @param {node} nTr the table row to check
		 *  @returns {boolean} true if the row is currently open, false otherwise
		 *  @dtopt API
		 *  @deprecated Since v1.10
		 *
		 *  @example
		 *    $(document).ready(function() {
		 *      var oTable;
		 *
		 *      // 'open' an information row when a row is clicked on
		 *      $('#example tbody tr').click( function () {
		 *        if ( oTable.fnIsOpen(this) ) {
		 *          oTable.fnClose( this );
		 *        } else {
		 *          oTable.fnOpen( this, "Temporary row opened", "info_row" );
		 *        }
		 *      } );
		 *
		 *      oTable = $('#example').dataTable();
		 *    } );
		 */
		this.fnIsOpen = function( nTr )
		{
			return this.api( true ).row( nTr ).child.isShown();
		};
		
		
		/**
		 * This function will place a new row directly after a row which is currently
		 * on display on the page, with the HTML contents that is passed into the
		 * function. This can be used, for example, to ask for confirmation that a
		 * particular record should be deleted.
		 *  @param {node} nTr The table row to 'open'
		 *  @param {string|node|jQuery} mHtml The HTML to put into the row
		 *  @param {string} sClass Class to give the new TD cell
		 *  @returns {node} The row opened. Note that if the table row passed in as the
		 *    first parameter, is not found in the table, this method will silently
		 *    return.
		 *  @dtopt API
		 *  @deprecated Since v1.10
		 *
		 *  @example
		 *    $(document).ready(function() {
		 *      var oTable;
		 *
		 *      // 'open' an information row when a row is clicked on
		 *      $('#example tbody tr').click( function () {
		 *        if ( oTable.fnIsOpen(this) ) {
		 *          oTable.fnClose( this );
		 *        } else {
		 *          oTable.fnOpen( this, "Temporary row opened", "info_row" );
		 *        }
		 *      } );
		 *
		 *      oTable = $('#example').dataTable();
		 *    } );
		 */
		this.fnOpen = function( nTr, mHtml, sClass )
		{
			return this.api( true )
				.row( nTr )
				.child( mHtml, sClass )
				.show()
				.child()[0];
		};
		
		
		/**
		 * Change the pagination - provides the internal logic for pagination in a simple API
		 * function. With this function you can have a DataTables table go to the next,
		 * previous, first or last pages.
		 *  @param {string|int} mAction Paging action to take: "first", "previous", "next" or "last"
		 *    or page number to jump to (integer), note that page 0 is the first page.
		 *  @param {bool} [bRedraw=true] Redraw the table or not
		 *  @dtopt API
		 *  @deprecated Since v1.10
		 *
		 *  @example
		 *    $(document).ready(function() {
		 *      var oTable = $('#example').dataTable();
		 *      oTable.fnPageChange( 'next' );
		 *    } );
		 */
		this.fnPageChange = function ( mAction, bRedraw )
		{
			var api = this.api( true ).page( mAction );
		
			if ( bRedraw === undefined || bRedraw ) {
				api.draw(false);
			}
		};
		
		
		/**
		 * Show a particular column
		 *  @param {int} iCol The column whose display should be changed
		 *  @param {bool} bShow Show (true) or hide (false) the column
		 *  @param {bool} [bRedraw=true] Redraw the table or not
		 *  @dtopt API
		 *  @deprecated Since v1.10
		 *
		 *  @example
		 *    $(document).ready(function() {
		 *      var oTable = $('#example').dataTable();
		 *
		 *      // Hide the second column after initialisation
		 *      oTable.fnSetColumnVis( 1, false );
		 *    } );
		 */
		this.fnSetColumnVis = function ( iCol, bShow, bRedraw )
		{
			var api = this.api( true ).column( iCol ).visible( bShow );
		
			if ( bRedraw === undefined || bRedraw ) {
				api.columns.adjust().draw();
			}
		};
		
		
		/**
		 * Get the settings for a particular table for external manipulation
		 *  @returns {object} DataTables settings object. See
		 *    {@link DataTable.models.oSettings}
		 *  @dtopt API
		 *  @deprecated Since v1.10
		 *
		 *  @example
		 *    $(document).ready(function() {
		 *      var oTable = $('#example').dataTable();
		 *      var oSettings = oTable.fnSettings();
		 *
		 *      // Show an example parameter from the settings
		 *      alert( oSettings._iDisplayStart );
		 *    } );
		 */
		this.fnSettings = function()
		{
			return _fnSettingsFromNode( this[_ext.iApiIndex] );
		};
		
		
		/**
		 * Sort the table by a particular column
		 *  @param {int} iCol the data index to sort on. Note that this will not match the
		 *    'display index' if you have hidden data entries
		 *  @dtopt API
		 *  @deprecated Since v1.10
		 *
		 *  @example
		 *    $(document).ready(function() {
		 *      var oTable = $('#example').dataTable();
		 *
		 *      // Sort immediately with columns 0 and 1
		 *      oTable.fnSort( [ [0,'asc'], [1,'asc'] ] );
		 *    } );
		 */
		this.fnSort = function( aaSort )
		{
			this.api( true ).order( aaSort ).draw();
		};
		
		
		/**
		 * Attach a sort listener to an element for a given column
		 *  @param {node} nNode the element to attach the sort listener to
		 *  @param {int} iColumn the column that a click on this node will sort on
		 *  @param {function} [fnCallback] callback function when sort is run
		 *  @dtopt API
		 *  @deprecated Since v1.10
		 *
		 *  @example
		 *    $(document).ready(function() {
		 *      var oTable = $('#example').dataTable();
		 *
		 *      // Sort on column 1, when 'sorter' is clicked on
		 *      oTable.fnSortListener( document.getElementById('sorter'), 1 );
		 *    } );
		 */
		this.fnSortListener = function( nNode, iColumn, fnCallback )
		{
			this.api( true ).order.listener( nNode, iColumn, fnCallback );
		};
		
		
		/**
		 * Update a table cell or row - this method will accept either a single value to
		 * update the cell with, an array of values with one element for each column or
		 * an object in the same format as the original data source. The function is
		 * self-referencing in order to make the multi column updates easier.
		 *  @param {object|array|string} mData Data to update the cell/row with
		 *  @param {node|int} mRow TR element you want to update or the aoData index
		 *  @param {int} [iColumn] The column to update, give as null or undefined to
		 *    update a whole row.
		 *  @param {bool} [bRedraw=true] Redraw the table or not
		 *  @param {bool} [bAction=true] Perform pre-draw actions or not
		 *  @returns {int} 0 on success, 1 on error
		 *  @dtopt API
		 *  @deprecated Since v1.10
		 *
		 *  @example
		 *    $(document).ready(function() {
		 *      var oTable = $('#example').dataTable();
		 *      oTable.fnUpdate( 'Example update', 0, 0 ); // Single cell
		 *      oTable.fnUpdate( ['a', 'b', 'c', 'd', 'e'], $('tbody tr')[0] ); // Row
		 *    } );
		 */
		this.fnUpdate = function( mData, mRow, iColumn, bRedraw, bAction )
		{
			var api = this.api( true );
		
			if ( iColumn === undefined || iColumn === null ) {
				api.row( mRow ).data( mData );
			}
			else {
				api.cell( mRow, iColumn ).data( mData );
			}
		
			if ( bAction === undefined || bAction ) {
				api.columns.adjust();
			}
		
			if ( bRedraw === undefined || bRedraw ) {
				api.draw();
			}
			return 0;
		};
		
		
		/**
		 * Provide a common method for plug-ins to check the version of DataTables being used, in order
		 * to ensure compatibility.
		 *  @param {string} sVersion Version string to check for, in the format "X.Y.Z". Note that the
		 *    formats "X" and "X.Y" are also acceptable.
		 *  @returns {boolean} true if this version of DataTables is greater or equal to the required
		 *    version, or false if this version of DataTales is not suitable
		 *  @method
		 *  @dtopt API
		 *  @deprecated Since v1.10
		 *
		 *  @example
		 *    $(document).ready(function() {
		 *      var oTable = $('#example').dataTable();
		 *      alert( oTable.fnVersionCheck( '1.9.0' ) );
		 *    } );
		 */
		this.fnVersionCheck = _ext.fnVersionCheck;
		
	
		var _that = this;
		var emptyInit = options === undefined;
		var len = this.length;
	
		if ( emptyInit ) {
			options = {};
		}
	
		this.oApi = this.internal = _ext.internal;
	
		// Extend with old style plug-in API methods
		for ( var fn in DataTable.ext.internal ) {
			if ( fn ) {
				this[fn] = _fnExternApiFunc(fn);
			}
		}
	
		this.each(function() {
			// For each initialisation we want to give it a clean initialisation
			// object that can be bashed around
			var o = {};
			var oInit = len > 1 ? // optimisation for single table case
				_fnExtend( o, options, true ) :
				options;
	
			/*global oInit,_that,emptyInit*/
			var i=0, iLen, j, jLen, k, kLen;
			var sId = this.getAttribute( 'id' );
			var bInitHandedOff = false;
			var defaults = DataTable.defaults;
			var $this = $(this);
			
			
			/* Sanity check */
			if ( this.nodeName.toLowerCase() != 'table' )
			{
				_fnLog( null, 0, 'Non-table node initialisation ('+this.nodeName+')', 2 );
				return;
			}
			
			/* Backwards compatibility for the defaults */
			_fnCompatOpts( defaults );
			_fnCompatCols( defaults.column );
			
			/* Convert the camel-case defaults to Hungarian */
			_fnCamelToHungarian( defaults, defaults, true );
			_fnCamelToHungarian( defaults.column, defaults.column, true );
			
			/* Setting up the initialisation object */
			_fnCamelToHungarian( defaults, $.extend( oInit, $this.data() ), true );
			
			
			
			/* Check to see if we are re-initialising a table */
			var allSettings = DataTable.settings;
			for ( i=0, iLen=allSettings.length ; i<iLen ; i++ )
			{
				var s = allSettings[i];
			
				/* Base check on table node */
				if (
					s.nTable == this ||
					(s.nTHead && s.nTHead.parentNode == this) ||
					(s.nTFoot && s.nTFoot.parentNode == this)
				) {
					var bRetrieve = oInit.bRetrieve !== undefined ? oInit.bRetrieve : defaults.bRetrieve;
					var bDestroy = oInit.bDestroy !== undefined ? oInit.bDestroy : defaults.bDestroy;
			
					if ( emptyInit || bRetrieve )
					{
						return s.oInstance;
					}
					else if ( bDestroy )
					{
						s.oInstance.fnDestroy();
						break;
					}
					else
					{
						_fnLog( s, 0, 'Cannot reinitialise DataTable', 3 );
						return;
					}
				}
			
				/* If the element we are initialising has the same ID as a table which was previously
				 * initialised, but the table nodes don't match (from before) then we destroy the old
				 * instance by simply deleting it. This is under the assumption that the table has been
				 * destroyed by other methods. Anyone using non-id selectors will need to do this manually
				 */
				if ( s.sTableId == this.id )
				{
					allSettings.splice( i, 1 );
					break;
				}
			}
			
			/* Ensure the table has an ID - required for accessibility */
			if ( sId === null || sId === "" )
			{
				sId = "DataTables_Table_"+(DataTable.ext._unique++);
				this.id = sId;
			}
			
			/* Create the settings object for this table and set some of the default parameters */
			var oSettings = $.extend( true, {}, DataTable.models.oSettings, {
				"sDestroyWidth": $this[0].style.width,
				"sInstance":     sId,
				"sTableId":      sId
			} );
			oSettings.nTable = this;
			oSettings.oApi   = _that.internal;
			oSettings.oInit  = oInit;
			
			allSettings.push( oSettings );
			
			// Need to add the instance after the instance after the settings object has been added
			// to the settings array, so we can self reference the table instance if more than one
			oSettings.oInstance = (_that.length===1) ? _that : $this.dataTable();
			
			// Backwards compatibility, before we apply all the defaults
			_fnCompatOpts( oInit );
			_fnLanguageCompat( oInit.oLanguage );
			
			// If the length menu is given, but the init display length is not, use the length menu
			if ( oInit.aLengthMenu && ! oInit.iDisplayLength )
			{
				oInit.iDisplayLength = Array.isArray( oInit.aLengthMenu[0] ) ?
					oInit.aLengthMenu[0][0] : oInit.aLengthMenu[0];
			}
			
			// Apply the defaults and init options to make a single init object will all
			// options defined from defaults and instance options.
			oInit = _fnExtend( $.extend( true, {}, defaults ), oInit );
			
			
			// Map the initialisation options onto the settings object
			_fnMap( oSettings.oFeatures, oInit, [
				"bPaginate",
				"bLengthChange",
				"bFilter",
				"bSort",
				"bSortMulti",
				"bInfo",
				"bProcessing",
				"bAutoWidth",
				"bSortClasses",
				"bServerSide",
				"bDeferRender"
			] );
			_fnMap( oSettings, oInit, [
				"asStripeClasses",
				"ajax",
				"fnServerData",
				"fnFormatNumber",
				"sServerMethod",
				"aaSorting",
				"aaSortingFixed",
				"aLengthMenu",
				"sPaginationType",
				"sAjaxSource",
				"sAjaxDataProp",
				"iStateDuration",
				"sDom",
				"bSortCellsTop",
				"iTabIndex",
				"fnStateLoadCallback",
				"fnStateSaveCallback",
				"renderer",
				"searchDelay",
				"rowId",
				[ "iCookieDuration", "iStateDuration" ], // backwards compat
				[ "oSearch", "oPreviousSearch" ],
				[ "aoSearchCols", "aoPreSearchCols" ],
				[ "iDisplayLength", "_iDisplayLength" ]
			] );
			_fnMap( oSettings.oScroll, oInit, [
				[ "sScrollX", "sX" ],
				[ "sScrollXInner", "sXInner" ],
				[ "sScrollY", "sY" ],
				[ "bScrollCollapse", "bCollapse" ]
			] );
			_fnMap( oSettings.oLanguage, oInit, "fnInfoCallback" );
			
			/* Callback functions which are array driven */
			_fnCallbackReg( oSettings, 'aoDrawCallback',       oInit.fnDrawCallback,      'user' );
			_fnCallbackReg( oSettings, 'aoServerParams',       oInit.fnServerParams,      'user' );
			_fnCallbackReg( oSettings, 'aoStateSaveParams',    oInit.fnStateSaveParams,   'user' );
			_fnCallbackReg( oSettings, 'aoStateLoadParams',    oInit.fnStateLoadParams,   'user' );
			_fnCallbackReg( oSettings, 'aoStateLoaded',        oInit.fnStateLoaded,       'user' );
			_fnCallbackReg( oSettings, 'aoRowCallback',        oInit.fnRowCallback,       'user' );
			_fnCallbackReg( oSettings, 'aoRowCreatedCallback', oInit.fnCreatedRow,        'user' );
			_fnCallbackReg( oSettings, 'aoHeaderCallback',     oInit.fnHeaderCallback,    'user' );
			_fnCallbackReg( oSettings, 'aoFooterCallback',     oInit.fnFooterCallback,    'user' );
			_fnCallbackReg( oSettings, 'aoInitComplete',       oInit.fnInitComplete,      'user' );
			_fnCallbackReg( oSettings, 'aoPreDrawCallback',    oInit.fnPreDrawCallback,   'user' );
			
			oSettings.rowIdFn = _fnGetObjectDataFn( oInit.rowId );
			
			/* Browser support detection */
			_fnBrowserDetect( oSettings );
			
			var oClasses = oSettings.oClasses;
			
			$.extend( oClasses, DataTable.ext.classes, oInit.oClasses );
			$this.addClass( oClasses.sTable );
			
			
			if ( oSettings.iInitDisplayStart === undefined )
			{
				/* Display start point, taking into account the save saving */
				oSettings.iInitDisplayStart = oInit.iDisplayStart;
				oSettings._iDisplayStart = oInit.iDisplayStart;
			}
			
			if ( oInit.iDeferLoading !== null )
			{
				oSettings.bDeferLoading = true;
				var tmp = Array.isArray( oInit.iDeferLoading );
				oSettings._iRecordsDisplay = tmp ? oInit.iDeferLoading[0] : oInit.iDeferLoading;
				oSettings._iRecordsTotal = tmp ? oInit.iDeferLoading[1] : oInit.iDeferLoading;
			}
			
			/* Language definitions */
			var oLanguage = oSettings.oLanguage;
			$.extend( true, oLanguage, oInit.oLanguage );
			
			if ( oLanguage.sUrl )
			{
				/* Get the language definitions from a file - because this Ajax call makes the language
				 * get async to the remainder of this function we use bInitHandedOff to indicate that
				 * _fnInitialise will be fired by the returned Ajax handler, rather than the constructor
				 */
				$.ajax( {
					dataType: 'json',
					url: oLanguage.sUrl,
					success: function ( json ) {
						_fnCamelToHungarian( defaults.oLanguage, json );
						_fnLanguageCompat( json );
						$.extend( true, oLanguage, json, oSettings.oInit.oLanguage );
			
						_fnCallbackFire( oSettings, null, 'i18n', [oSettings]);
						_fnInitialise( oSettings );
					},
					error: function () {
						// Error occurred loading language file, continue on as best we can
						_fnInitialise( oSettings );
					}
				} );
				bInitHandedOff = true;
			}
			else {
				_fnCallbackFire( oSettings, null, 'i18n', [oSettings]);
			}
			
			/*
			 * Stripes
			 */
			if ( oInit.asStripeClasses === null )
			{
				oSettings.asStripeClasses =[
					oClasses.sStripeOdd,
					oClasses.sStripeEven
				];
			}
			
			/* Remove row stripe classes if they are already on the table row */
			var stripeClasses = oSettings.asStripeClasses;
			var rowOne = $this.children('tbody').find('tr').eq(0);
			if ( $.inArray( true, $.map( stripeClasses, function(el, i) {
				return rowOne.hasClass(el);
			} ) ) !== -1 ) {
				$('tbody tr', this).removeClass( stripeClasses.join(' ') );
				oSettings.asDestroyStripes = stripeClasses.slice();
			}
			
			/*
			 * Columns
			 * See if we should load columns automatically or use defined ones
			 */
			var anThs = [];
			var aoColumnsInit;
			var nThead = this.getElementsByTagName('thead');
			if ( nThead.length !== 0 )
			{
				_fnDetectHeader( oSettings.aoHeader, nThead[0] );
				anThs = _fnGetUniqueThs( oSettings );
			}
			
			/* If not given a column array, generate one with nulls */
			if ( oInit.aoColumns === null )
			{
				aoColumnsInit = [];
				for ( i=0, iLen=anThs.length ; i<iLen ; i++ )
				{
					aoColumnsInit.push( null );
				}
			}
			else
			{
				aoColumnsInit = oInit.aoColumns;
			}
			
			/* Add the columns */
			for ( i=0, iLen=aoColumnsInit.length ; i<iLen ; i++ )
			{
				_fnAddColumn( oSettings, anThs ? anThs[i] : null );
			}
			
			/* Apply the column definitions */
			_fnApplyColumnDefs( oSettings, oInit.aoColumnDefs, aoColumnsInit, function (iCol, oDef) {
				_fnColumnOptions( oSettings, iCol, oDef );
			} );
			
			/* HTML5 attribute detection - build an mData object automatically if the
			 * attributes are found
			 */
			if ( rowOne.length ) {
				var a = function ( cell, name ) {
					return cell.getAttribute( 'data-'+name ) !== null ? name : null;
				};
			
				$( rowOne[0] ).children('th, td').each( function (i, cell) {
					var col = oSettings.aoColumns[i];
			
					if (! col) {
						_fnLog( oSettings, 0, 'Incorrect column count', 18 );
					}
			
					if ( col.mData === i ) {
						var sort = a( cell, 'sort' ) || a( cell, 'order' );
						var filter = a( cell, 'filter' ) || a( cell, 'search' );
			
						if ( sort !== null || filter !== null ) {
							col.mData = {
								_:      i+'.display',
								sort:   sort !== null   ? i+'.@data-'+sort   : undefined,
								type:   sort !== null   ? i+'.@data-'+sort   : undefined,
								filter: filter !== null ? i+'.@data-'+filter : undefined
							};
							col._isArrayHost = true;
			
							_fnColumnOptions( oSettings, i );
						}
					}
				} );
			}
			
			var features = oSettings.oFeatures;
			var loadedInit = function () {
				/*
				 * Sorting
				 * @todo For modularisation (1.11) this needs to do into a sort start up handler
				 */
			
				// If aaSorting is not defined, then we use the first indicator in asSorting
				// in case that has been altered, so the default sort reflects that option
				if ( oInit.aaSorting === undefined ) {
					var sorting = oSettings.aaSorting;
					for ( i=0, iLen=sorting.length ; i<iLen ; i++ ) {
						sorting[i][1] = oSettings.aoColumns[ i ].asSorting[0];
					}
				}
			
				/* Do a first pass on the sorting classes (allows any size changes to be taken into
				 * account, and also will apply sorting disabled classes if disabled
				 */
				_fnSortingClasses( oSettings );
			
				if ( features.bSort ) {
					_fnCallbackReg( oSettings, 'aoDrawCallback', function () {
						if ( oSettings.bSorted ) {
							var aSort = _fnSortFlatten( oSettings );
							var sortedColumns = {};
			
							$.each( aSort, function (i, val) {
								sortedColumns[ val.src ] = val.dir;
							} );
			
							_fnCallbackFire( oSettings, null, 'order', [oSettings, aSort, sortedColumns] );
							_fnSortAria( oSettings );
						}
					} );
				}
			
				_fnCallbackReg( oSettings, 'aoDrawCallback', function () {
					if ( oSettings.bSorted || _fnDataSource( oSettings ) === 'ssp' || features.bDeferRender ) {
						_fnSortingClasses( oSettings );
					}
				}, 'sc' );
			
			
				/*
				 * Final init
				 * Cache the header, body and footer as required, creating them if needed
				 */
			
				// Work around for Webkit bug 83867 - store the caption-side before removing from doc
				var captions = $this.children('caption').each( function () {
					this._captionSide = $(this).css('caption-side');
				} );
			
				var thead = $this.children('thead');
				if ( thead.length === 0 ) {
					thead = $('<thead/>').appendTo($this);
				}
				oSettings.nTHead = thead[0];
			
				var tbody = $this.children('tbody');
				if ( tbody.length === 0 ) {
					tbody = $('<tbody/>').insertAfter(thead);
				}
				oSettings.nTBody = tbody[0];
			
				var tfoot = $this.children('tfoot');
				if ( tfoot.length === 0 && captions.length > 0 && (oSettings.oScroll.sX !== "" || oSettings.oScroll.sY !== "") ) {
					// If we are a scrolling table, and no footer has been given, then we need to create
					// a tfoot element for the caption element to be appended to
					tfoot = $('<tfoot/>').appendTo($this);
				}
			
				if ( tfoot.length === 0 || tfoot.children().length === 0 ) {
					$this.addClass( oClasses.sNoFooter );
				}
				else if ( tfoot.length > 0 ) {
					oSettings.nTFoot = tfoot[0];
					_fnDetectHeader( oSettings.aoFooter, oSettings.nTFoot );
				}
			
				/* Check if there is data passing into the constructor */
				if ( oInit.aaData ) {
					for ( i=0 ; i<oInit.aaData.length ; i++ ) {
						_fnAddData( oSettings, oInit.aaData[ i ] );
					}
				}
				else if ( oSettings.bDeferLoading || _fnDataSource( oSettings ) == 'dom' ) {
					/* Grab the data from the page - only do this when deferred loading or no Ajax
					 * source since there is no point in reading the DOM data if we are then going
					 * to replace it with Ajax data
					 */
					_fnAddTr( oSettings, $(oSettings.nTBody).children('tr') );
				}
			
				/* Copy the data index array */
				oSettings.aiDisplay = oSettings.aiDisplayMaster.slice();
			
				/* Initialisation complete - table can be drawn */
				oSettings.bInitialised = true;
			
				/* Check if we need to initialise the table (it might not have been handed off to the
				 * language processor)
				 */
				if ( bInitHandedOff === false ) {
					_fnInitialise( oSettings );
				}
			};
			
			/* Must be done after everything which can be overridden by the state saving! */
			_fnCallbackReg( oSettings, 'aoDrawCallback', _fnSaveState, 'state_save' );
			
			if ( oInit.bStateSave )
			{
				features.bStateSave = true;
				_fnLoadState( oSettings, oInit, loadedInit );
			}
			else {
				loadedInit();
			}
			
		} );
		_that = null;
		return this;
	};
	
	
	/*
	 * It is useful to have variables which are scoped locally so only the
	 * DataTables functions can access them and they don't leak into global space.
	 * At the same time these functions are often useful over multiple files in the
	 * core and API, so we list, or at least document, all variables which are used
	 * by DataTables as private variables here. This also ensures that there is no
	 * clashing of variable names and that they can easily referenced for reuse.
	 */
	
	
	// Defined else where
	//  _selector_run
	//  _selector_opts
	//  _selector_first
	//  _selector_row_indexes
	
	var _ext; // DataTable.ext
	var _Api; // DataTable.Api
	var _api_register; // DataTable.Api.register
	var _api_registerPlural; // DataTable.Api.registerPlural
	
	var _re_dic = {};
	var _re_new_lines = /[\r\n\u2028]/g;
	var _re_html = /<.*?>/g;
	
	// This is not strict ISO8601 - Date.parse() is quite lax, although
	// implementations differ between browsers.
	var _re_date = /^\d{2,4}[\.\/\-]\d{1,2}[\.\/\-]\d{1,2}([T ]{1}\d{1,2}[:\.]\d{2}([\.:]\d{2})?)?$/;
	
	// Escape regular expression special characters
	var _re_escape_regex = new RegExp( '(\\' + [ '/', '.', '*', '+', '?', '|', '(', ')', '[', ']', '{', '}', '\\', '$', '^', '-' ].join('|\\') + ')', 'g' );
	
	// http://en.wikipedia.org/wiki/Foreign_exchange_market
	// - \u20BD - Russian ruble.
	// - \u20a9 - South Korean Won
	// - \u20BA - Turkish Lira
	// - \u20B9 - Indian Rupee
	// - R - Brazil (R$) and South Africa
	// - fr - Swiss Franc
	// - kr - Swedish krona, Norwegian krone and Danish krone
	// - \u2009 is thin space and \u202F is narrow no-break space, both used in many
	// - Ƀ - Bitcoin
	// - Ξ - Ethereum
	//   standards as thousands separators.
	var _re_formatted_numeric = /['\u00A0,$£€¥%\u2009\u202F\u20BD\u20a9\u20BArfkɃΞ]/gi;
	
	
	var _empty = function ( d ) {
		return !d || d === true || d === '-' ? true : false;
	};
	
	
	var _intVal = function ( s ) {
		var integer = parseInt( s, 10 );
		return !isNaN(integer) && isFinite(s) ? integer : null;
	};
	
	// Convert from a formatted number with characters other than `.` as the
	// decimal place, to a Javascript number
	var _numToDecimal = function ( num, decimalPoint ) {
		// Cache created regular expressions for speed as this function is called often
		if ( ! _re_dic[ decimalPoint ] ) {
			_re_dic[ decimalPoint ] = new RegExp( _fnEscapeRegex( decimalPoint ), 'g' );
		}
		return typeof num === 'string' && decimalPoint !== '.' ?
			num.replace( /\./g, '' ).replace( _re_dic[ decimalPoint ], '.' ) :
			num;
	};
	
	
	var _isNumber = function ( d, decimalPoint, formatted ) {
		let type = typeof d;
		var strType = type === 'string';
	
		if ( type === 'number' || type === 'bigint') {
			return true;
		}
	
		// If empty return immediately so there must be a number if it is a
		// formatted string (this stops the string "k", or "kr", etc being detected
		// as a formatted number for currency
		if ( _empty( d ) ) {
			return true;
		}
	
		if ( decimalPoint && strType ) {
			d = _numToDecimal( d, decimalPoint );
		}
	
		if ( formatted && strType ) {
			d = d.replace( _re_formatted_numeric, '' );
		}
	
		return !isNaN( parseFloat(d) ) && isFinite( d );
	};
	
	
	// A string without HTML in it can be considered to be HTML still
	var _isHtml = function ( d ) {
		return _empty( d ) || typeof d === 'string';
	};
	
	
	var _htmlNumeric = function ( d, decimalPoint, formatted ) {
		if ( _empty( d ) ) {
			return true;
		}
	
		var html = _isHtml( d );
		return ! html ?
			null :
			_isNumber( _stripHtml( d ), decimalPoint, formatted ) ?
				true :
				null;
	};
	
	
	var _pluck = function ( a, prop, prop2 ) {
		var out = [];
		var i=0, ien=a.length;
	
		// Could have the test in the loop for slightly smaller code, but speed
		// is essential here
		if ( prop2 !== undefined ) {
			for ( ; i<ien ; i++ ) {
				if ( a[i] && a[i][ prop ] ) {
					out.push( a[i][ prop ][ prop2 ] );
				}
			}
		}
		else {
			for ( ; i<ien ; i++ ) {
				if ( a[i] ) {
					out.push( a[i][ prop ] );
				}
			}
		}
	
		return out;
	};
	
	
	// Basically the same as _pluck, but rather than looping over `a` we use `order`
	// as the indexes to pick from `a`
	var _pluck_order = function ( a, order, prop, prop2 )
	{
		var out = [];
		var i=0, ien=order.length;
	
		// Could have the test in the loop for slightly smaller code, but speed
		// is essential here
		if ( prop2 !== undefined ) {
			for ( ; i<ien ; i++ ) {
				if ( a[ order[i] ][ prop ] ) {
					out.push( a[ order[i] ][ prop ][ prop2 ] );
				}
			}
		}
		else {
			for ( ; i<ien ; i++ ) {
				out.push( a[ order[i] ][ prop ] );
			}
		}
	
		return out;
	};
	
	
	var _range = function ( len, start )
	{
		var out = [];
		var end;
	
		if ( start === undefined ) {
			start = 0;
			end = len;
		}
		else {
			end = start;
			start = len;
		}
	
		for ( var i=start ; i<end ; i++ ) {
			out.push( i );
		}
	
		return out;
	};
	
	
	var _removeEmpty = function ( a )
	{
		var out = [];
	
		for ( var i=0, ien=a.length ; i<ien ; i++ ) {
			if ( a[i] ) { // careful - will remove all falsy values!
				out.push( a[i] );
			}
		}
	
		return out;
	};
	
	
	var _stripHtml = function ( d ) {
		return d.replace( _re_html, '' );
	};
	
	
	/**
	 * Determine if all values in the array are unique. This means we can short
	 * cut the _unique method at the cost of a single loop. A sorted array is used
	 * to easily check the values.
	 *
	 * @param  {array} src Source array
	 * @return {boolean} true if all unique, false otherwise
	 * @ignore
	 */
	var _areAllUnique = function ( src ) {
		if ( src.length < 2 ) {
			return true;
		}
	
		var sorted = src.slice().sort();
		var last = sorted[0];
	
		for ( var i=1, ien=sorted.length ; i<ien ; i++ ) {
			if ( sorted[i] === last ) {
				return false;
			}
	
			last = sorted[i];
		}
	
		return true;
	};
	
	
	/**
	 * Find the unique elements in a source array.
	 *
	 * @param  {array} src Source array
	 * @return {array} Array of unique items
	 * @ignore
	 */
	var _unique = function ( src )
	{
		if ( _areAllUnique( src ) ) {
			return src.slice();
		}
	
		// A faster unique method is to use object keys to identify used values,
		// but this doesn't work with arrays or objects, which we must also
		// consider. See jsperf.com/compare-array-unique-versions/4 for more
		// information.
		var
			out = [],
			val,
			i, ien=src.length,
			j, k=0;
	
		again: for ( i=0 ; i<ien ; i++ ) {
			val = src[i];
	
			for ( j=0 ; j<k ; j++ ) {
				if ( out[j] === val ) {
					continue again;
				}
			}
	
			out.push( val );
			k++;
		}
	
		return out;
	};
	
	// Surprisingly this is faster than [].concat.apply
	// https://jsperf.com/flatten-an-array-loop-vs-reduce/2
	var _flatten = function (out, val) {
		if (Array.isArray(val)) {
			for (var i=0 ; i<val.length ; i++) {
				_flatten(out, val[i]);
			}
		}
		else {
			out.push(val);
		}
	  
		return out;
	}
	
	var _includes = function (search, start) {
		if (start === undefined) {
			start = 0;
		}
	
		return this.indexOf(search, start) !== -1;	
	};
	
	// Array.isArray polyfill.
	// https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Array/isArray
	if (! Array.isArray) {
	    Array.isArray = function(arg) {
	        return Object.prototype.toString.call(arg) === '[object Array]';
	    };
	}
	
	if (! Array.prototype.includes) {
		Array.prototype.includes = _includes;
	}
	
	// .trim() polyfill
	// https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/String/trim
	if (!String.prototype.trim) {
	  String.prototype.trim = function () {
	    return this.replace(/^[\s\uFEFF\xA0]+|[\s\uFEFF\xA0]+$/g, '');
	  };
	}
	
	if (! String.prototype.includes) {
		String.prototype.includes = _includes;
	}
	
	/**
	 * DataTables utility methods
	 * 
	 * This namespace provides helper methods that DataTables uses internally to
	 * create a DataTable, but which are not exclusively used only for DataTables.
	 * These methods can be used by extension authors to save the duplication of
	 * code.
	 *
	 *  @namespace
	 */
	DataTable.util = {
		/**
		 * Throttle the calls to a function. Arguments and context are maintained
		 * for the throttled function.
		 *
		 * @param {function} fn Function to be called
		 * @param {integer} freq Call frequency in mS
		 * @return {function} Wrapped function
		 */
		throttle: function ( fn, freq ) {
			var
				frequency = freq !== undefined ? freq : 200,
				last,
				timer;
	
			return function () {
				var
					that = this,
					now  = +new Date(),
					args = arguments;
	
				if ( last && now < last + frequency ) {
					clearTimeout( timer );
	
					timer = setTimeout( function () {
						last = undefined;
						fn.apply( that, args );
					}, frequency );
				}
				else {
					last = now;
					fn.apply( that, args );
				}
			};
		},
	
	
		/**
		 * Escape a string such that it can be used in a regular expression
		 *
		 *  @param {string} val string to escape
		 *  @returns {string} escaped string
		 */
		escapeRegex: function ( val ) {
			return val.replace( _re_escape_regex, '\\$1' );
		},
	
		/**
		 * Create a function that will write to a nested object or array
		 * @param {*} source JSON notation string
		 * @returns Write function
		 */
		set: function ( source ) {
			if ( $.isPlainObject( source ) ) {
				/* Unlike get, only the underscore (global) option is used for for
				 * setting data since we don't know the type here. This is why an object
				 * option is not documented for `mData` (which is read/write), but it is
				 * for `mRender` which is read only.
				 */
				return DataTable.util.set( source._ );
			}
			else if ( source === null ) {
				// Nothing to do when the data source is null
				return function () {};
			}
			else if ( typeof source === 'function' ) {
				return function (data, val, meta) {
					source( data, 'set', val, meta );
				};
			}
			else if ( typeof source === 'string' && (source.indexOf('.') !== -1 ||
					  source.indexOf('[') !== -1 || source.indexOf('(') !== -1) )
			{
				// Like the get, we need to get data from a nested object
				var setData = function (data, val, src) {
					var a = _fnSplitObjNotation( src ), b;
					var aLast = a[a.length-1];
					var arrayNotation, funcNotation, o, innerSrc;
		
					for ( var i=0, iLen=a.length-1 ; i<iLen ; i++ ) {
						// Protect against prototype pollution
						if (a[i] === '__proto__' || a[i] === 'constructor') {
							throw new Error('Cannot set prototype values');
						}
		
						// Check if we are dealing with an array notation request
						arrayNotation = a[i].match(__reArray);
						funcNotation = a[i].match(__reFn);
		
						if ( arrayNotation ) {
							a[i] = a[i].replace(__reArray, '');
							data[ a[i] ] = [];
		
							// Get the remainder of the nested object to set so we can recurse
							b = a.slice();
							b.splice( 0, i+1 );
							innerSrc = b.join('.');
		
							// Traverse each entry in the array setting the properties requested
							if ( Array.isArray( val ) ) {
								for ( var j=0, jLen=val.length ; j<jLen ; j++ ) {
									o = {};
									setData( o, val[j], innerSrc );
									data[ a[i] ].push( o );
								}
							}
							else {
								// We've been asked to save data to an array, but it
								// isn't array data to be saved. Best that can be done
								// is to just save the value.
								data[ a[i] ] = val;
							}
		
							// The inner call to setData has already traversed through the remainder
							// of the source and has set the data, thus we can exit here
							return;
						}
						else if ( funcNotation ) {
							// Function call
							a[i] = a[i].replace(__reFn, '');
							data = data[ a[i] ]( val );
						}
		
						// If the nested object doesn't currently exist - since we are
						// trying to set the value - create it
						if ( data[ a[i] ] === null || data[ a[i] ] === undefined ) {
							data[ a[i] ] = {};
						}
						data = data[ a[i] ];
					}
		
					// Last item in the input - i.e, the actual set
					if ( aLast.match(__reFn ) ) {
						// Function call
						data = data[ aLast.replace(__reFn, '') ]( val );
					}
					else {
						// If array notation is used, we just want to strip it and use the property name
						// and assign the value. If it isn't used, then we get the result we want anyway
						data[ aLast.replace(__reArray, '') ] = val;
					}
				};
		
				return function (data, val) { // meta is also passed in, but not used
					return setData( data, val, source );
				};
			}
			else {
				// Array or flat object mapping
				return function (data, val) { // meta is also passed in, but not used
					data[source] = val;
				};
			}
		},
	
		/**
		 * Create a function that will read nested objects from arrays, based on JSON notation
		 * @param {*} source JSON notation string
		 * @returns Value read
		 */
		get: function ( source ) {
			if ( $.isPlainObject( source ) ) {
				// Build an object of get functions, and wrap them in a single call
				var o = {};
				$.each( source, function (key, val) {
					if ( val ) {
						o[key] = DataTable.util.get( val );
					}
				} );
		
				return function (data, type, row, meta) {
					var t = o[type] || o._;
					return t !== undefined ?
						t(data, type, row, meta) :
						data;
				};
			}
			else if ( source === null ) {
				// Give an empty string for rendering / sorting etc
				return function (data) { // type, row and meta also passed, but not used
					return data;
				};
			}
			else if ( typeof source === 'function' ) {
				return function (data, type, row, meta) {
					return source( data, type, row, meta );
				};
			}
			else if ( typeof source === 'string' && (source.indexOf('.') !== -1 ||
					  source.indexOf('[') !== -1 || source.indexOf('(') !== -1) )
			{
				/* If there is a . in the source string then the data source is in a
				 * nested object so we loop over the data for each level to get the next
				 * level down. On each loop we test for undefined, and if found immediately
				 * return. This allows entire objects to be missing and sDefaultContent to
				 * be used if defined, rather than throwing an error
				 */
				var fetchData = function (data, type, src) {
					var arrayNotation, funcNotation, out, innerSrc;
		
					if ( src !== "" ) {
						var a = _fnSplitObjNotation( src );
		
						for ( var i=0, iLen=a.length ; i<iLen ; i++ ) {
							// Check if we are dealing with special notation
							arrayNotation = a[i].match(__reArray);
							funcNotation = a[i].match(__reFn);
		
							if ( arrayNotation ) {
								// Array notation
								a[i] = a[i].replace(__reArray, '');
		
								// Condition allows simply [] to be passed in
								if ( a[i] !== "" ) {
									data = data[ a[i] ];
								}
								out = [];
		
								// Get the remainder of the nested object to get
								a.splice( 0, i+1 );
								innerSrc = a.join('.');
		
								// Traverse each entry in the array getting the properties requested
								if ( Array.isArray( data ) ) {
									for ( var j=0, jLen=data.length ; j<jLen ; j++ ) {
										out.push( fetchData( data[j], type, innerSrc ) );
									}
								}
		
								// If a string is given in between the array notation indicators, that
								// is used to join the strings together, otherwise an array is returned
								var join = arrayNotation[0].substring(1, arrayNotation[0].length-1);
								data = (join==="") ? out : out.join(join);
		
								// The inner call to fetchData has already traversed through the remainder
								// of the source requested, so we exit from the loop
								break;
							}
							else if ( funcNotation ) {
								// Function call
								a[i] = a[i].replace(__reFn, '');
								data = data[ a[i] ]();
								continue;
							}
		
							if ( data === null || data[ a[i] ] === undefined ) {
								return undefined;
							}
	
							data = data[ a[i] ];
						}
					}
		
					return data;
				};
		
				return function (data, type) { // row and meta also passed, but not used
					return fetchData( data, type, source );
				};
			}
			else {
				// Array or flat object mapping
				return function (data, type) { // row and meta also passed, but not used
					return data[source];
				};
			}
		}
	};
	
	
	
	/**
	 * Create a mapping object that allows camel case parameters to be looked up
	 * for their Hungarian counterparts. The mapping is stored in a private
	 * parameter called `_hungarianMap` which can be accessed on the source object.
	 *  @param {object} o
	 *  @memberof DataTable#oApi
	 */
	function _fnHungarianMap ( o )
	{
		var
			hungarian = 'a aa ai ao as b fn i m o s ',
			match,
			newKey,
			map = {};
	
		$.each( o, function (key, val) {
			match = key.match(/^([^A-Z]+?)([A-Z])/);
	
			if ( match && hungarian.indexOf(match[1]+' ') !== -1 )
			{
				newKey = key.replace( match[0], match[2].toLowerCase() );
				map[ newKey ] = key;
	
				if ( match[1] === 'o' )
				{
					_fnHungarianMap( o[key] );
				}
			}
		} );
	
		o._hungarianMap = map;
	}
	
	
	/**
	 * Convert from camel case parameters to Hungarian, based on a Hungarian map
	 * created by _fnHungarianMap.
	 *  @param {object} src The model object which holds all parameters that can be
	 *    mapped.
	 *  @param {object} user The object to convert from camel case to Hungarian.
	 *  @param {boolean} force When set to `true`, properties which already have a
	 *    Hungarian value in the `user` object will be overwritten. Otherwise they
	 *    won't be.
	 *  @memberof DataTable#oApi
	 */
	function _fnCamelToHungarian ( src, user, force )
	{
		if ( ! src._hungarianMap ) {
			_fnHungarianMap( src );
		}
	
		var hungarianKey;
	
		$.each( user, function (key, val) {
			hungarianKey = src._hungarianMap[ key ];
	
			if ( hungarianKey !== undefined && (force || user[hungarianKey] === undefined) )
			{
				// For objects, we need to buzz down into the object to copy parameters
				if ( hungarianKey.charAt(0) === 'o' )
				{
					// Copy the camelCase options over to the hungarian
					if ( ! user[ hungarianKey ] ) {
						user[ hungarianKey ] = {};
					}
					$.extend( true, user[hungarianKey], user[key] );
	
					_fnCamelToHungarian( src[hungarianKey], user[hungarianKey], force );
				}
				else {
					user[hungarianKey] = user[ key ];
				}
			}
		} );
	}
	
	
	/**
	 * Language compatibility - when certain options are given, and others aren't, we
	 * need to duplicate the values over, in order to provide backwards compatibility
	 * with older language files.
	 *  @param {object} oSettings dataTables settings object
	 *  @memberof DataTable#oApi
	 */
	function _fnLanguageCompat( lang )
	{
		// Note the use of the Hungarian notation for the parameters in this method as
		// this is called after the mapping of camelCase to Hungarian
		var defaults = DataTable.defaults.oLanguage;
	
		// Default mapping
		var defaultDecimal = defaults.sDecimal;
		if ( defaultDecimal ) {
			_addNumericSort( defaultDecimal );
		}
	
		if ( lang ) {
			var zeroRecords = lang.sZeroRecords;
	
			// Backwards compatibility - if there is no sEmptyTable given, then use the same as
			// sZeroRecords - assuming that is given.
			if ( ! lang.sEmptyTable && zeroRecords &&
				defaults.sEmptyTable === "No data available in table" )
			{
				_fnMap( lang, lang, 'sZeroRecords', 'sEmptyTable' );
			}
	
			// Likewise with loading records
			if ( ! lang.sLoadingRecords && zeroRecords &&
				defaults.sLoadingRecords === "Loading..." )
			{
				_fnMap( lang, lang, 'sZeroRecords', 'sLoadingRecords' );
			}
	
			// Old parameter name of the thousands separator mapped onto the new
			if ( lang.sInfoThousands ) {
				lang.sThousands = lang.sInfoThousands;
			}
	
			var decimal = lang.sDecimal;
			if ( decimal && defaultDecimal !== decimal ) {
				_addNumericSort( decimal );
			}
		}
	}
	
	
	/**
	 * Map one parameter onto another
	 *  @param {object} o Object to map
	 *  @param {*} knew The new parameter name
	 *  @param {*} old The old parameter name
	 */
	var _fnCompatMap = function ( o, knew, old ) {
		if ( o[ knew ] !== undefined ) {
			o[ old ] = o[ knew ];
		}
	};
	
	
	/**
	 * Provide backwards compatibility for the main DT options. Note that the new
	 * options are mapped onto the old parameters, so this is an external interface
	 * change only.
	 *  @param {object} init Object to map
	 */
	function _fnCompatOpts ( init )
	{
		_fnCompatMap( init, 'ordering',      'bSort' );
		_fnCompatMap( init, 'orderMulti',    'bSortMulti' );
		_fnCompatMap( init, 'orderClasses',  'bSortClasses' );
		_fnCompatMap( init, 'orderCellsTop', 'bSortCellsTop' );
		_fnCompatMap( init, 'order',         'aaSorting' );
		_fnCompatMap( init, 'orderFixed',    'aaSortingFixed' );
		_fnCompatMap( init, 'paging',        'bPaginate' );
		_fnCompatMap( init, 'pagingType',    'sPaginationType' );
		_fnCompatMap( init, 'pageLength',    'iDisplayLength' );
		_fnCompatMap( init, 'searching',     'bFilter' );
	
		// Boolean initialisation of x-scrolling
		if ( typeof init.sScrollX === 'boolean' ) {
			init.sScrollX = init.sScrollX ? '100%' : '';
		}
		if ( typeof init.scrollX === 'boolean' ) {
			init.scrollX = init.scrollX ? '100%' : '';
		}
	
		// Column search objects are in an array, so it needs to be converted
		// element by element
		var searchCols = init.aoSearchCols;
	
		if ( searchCols ) {
			for ( var i=0, ien=searchCols.length ; i<ien ; i++ ) {
				if ( searchCols[i] ) {
					_fnCamelToHungarian( DataTable.models.oSearch, searchCols[i] );
				}
			}
		}
	}
	
	
	/**
	 * Provide backwards compatibility for column options. Note that the new options
	 * are mapped onto the old parameters, so this is an external interface change
	 * only.
	 *  @param {object} init Object to map
	 */
	function _fnCompatCols ( init )
	{
		_fnCompatMap( init, 'orderable',     'bSortable' );
		_fnCompatMap( init, 'orderData',     'aDataSort' );
		_fnCompatMap( init, 'orderSequence', 'asSorting' );
		_fnCompatMap( init, 'orderDataType', 'sortDataType' );
	
		// orderData can be given as an integer
		var dataSort = init.aDataSort;
		if ( typeof dataSort === 'number' && ! Array.isArray( dataSort ) ) {
			init.aDataSort = [ dataSort ];
		}
	}
	
	
	/**
	 * Browser feature detection for capabilities, quirks
	 *  @param {object} settings dataTables settings object
	 *  @memberof DataTable#oApi
	 */
	function _fnBrowserDetect( settings )
	{
		// We don't need to do this every time DataTables is constructed, the values
		// calculated are specific to the browser and OS configuration which we
		// don't expect to change between initialisations
		if ( ! DataTable.__browser ) {
			var browser = {};
			DataTable.__browser = browser;
	
			// Scrolling feature / quirks detection
			var n = $('<div/>')
				.css( {
					position: 'fixed',
					top: 0,
					left: $(window).scrollLeft()*-1, // allow for scrolling
					height: 1,
					width: 1,
					overflow: 'hidden'
				} )
				.append(
					$('<div/>')
						.css( {
							position: 'absolute',
							top: 1,
							left: 1,
							width: 100,
							overflow: 'scroll'
						} )
						.append(
							$('<div/>')
								.css( {
									width: '100%',
									height: 10
								} )
						)
				)
				.appendTo( 'body' );
	
			var outer = n.children();
			var inner = outer.children();
	
			// Numbers below, in order, are:
			// inner.offsetWidth, inner.clientWidth, outer.offsetWidth, outer.clientWidth
			//
			// IE6 XP:                           100 100 100  83
			// IE7 Vista:                        100 100 100  83
			// IE 8+ Windows:                     83  83 100  83
			// Evergreen Windows:                 83  83 100  83
			// Evergreen Mac with scrollbars:     85  85 100  85
			// Evergreen Mac without scrollbars: 100 100 100 100
	
			// Get scrollbar width
			browser.barWidth = outer[0].offsetWidth - outer[0].clientWidth;
	
			// IE6/7 will oversize a width 100% element inside a scrolling element, to
			// include the width of the scrollbar, while other browsers ensure the inner
			// element is contained without forcing scrolling
			browser.bScrollOversize = inner[0].offsetWidth === 100 && outer[0].clientWidth !== 100;
	
			// In rtl text layout, some browsers (most, but not all) will place the
			// scrollbar on the left, rather than the right.
			browser.bScrollbarLeft = Math.round( inner.offset().left ) !== 1;
	
			// IE8- don't provide height and width for getBoundingClientRect
			browser.bBounding = n[0].getBoundingClientRect().width ? true : false;
	
			n.remove();
		}
	
		$.extend( settings.oBrowser, DataTable.__browser );
		settings.oScroll.iBarWidth = DataTable.__browser.barWidth;
	}
	
	
	/**
	 * Array.prototype reduce[Right] method, used for browsers which don't support
	 * JS 1.6. Done this way to reduce code size, since we iterate either way
	 *  @param {object} settings dataTables settings object
	 *  @memberof DataTable#oApi
	 */
	function _fnReduce ( that, fn, init, start, end, inc )
	{
		var
			i = start,
			value,
			isSet = false;
	
		if ( init !== undefined ) {
			value = init;
			isSet = true;
		}
	
		while ( i !== end ) {
			if ( ! that.hasOwnProperty(i) ) {
				continue;
			}
	
			value = isSet ?
				fn( value, that[i], i, that ) :
				that[i];
	
			isSet = true;
			i += inc;
		}
	
		return value;
	}
	
	/**
	 * Add a column to the list used for the table with default values
	 *  @param {object} oSettings dataTables settings object
	 *  @param {node} nTh The th element for this column
	 *  @memberof DataTable#oApi
	 */
	function _fnAddColumn( oSettings, nTh )
	{
		// Add column to aoColumns array
		var oDefaults = DataTable.defaults.column;
		var iCol = oSettings.aoColumns.length;
		var oCol = $.extend( {}, DataTable.models.oColumn, oDefaults, {
			"nTh": nTh ? nTh : document.createElement('th'),
			"sTitle":    oDefaults.sTitle    ? oDefaults.sTitle    : nTh ? nTh.innerHTML : '',
			"aDataSort": oDefaults.aDataSort ? oDefaults.aDataSort : [iCol],
			"mData": oDefaults.mData ? oDefaults.mData : iCol,
			idx: iCol
		} );
		oSettings.aoColumns.push( oCol );
	
		// Add search object for column specific search. Note that the `searchCols[ iCol ]`
		// passed into extend can be undefined. This allows the user to give a default
		// with only some of the parameters defined, and also not give a default
		var searchCols = oSettings.aoPreSearchCols;
		searchCols[ iCol ] = $.extend( {}, DataTable.models.oSearch, searchCols[ iCol ] );
	
		// Use the default column options function to initialise classes etc
		_fnColumnOptions( oSettings, iCol, $(nTh).data() );
	}
	
	
	/**
	 * Apply options for a column
	 *  @param {object} oSettings dataTables settings object
	 *  @param {int} iCol column index to consider
	 *  @param {object} oOptions object with sType, bVisible and bSearchable etc
	 *  @memberof DataTable#oApi
	 */
	function _fnColumnOptions( oSettings, iCol, oOptions )
	{
		var oCol = oSettings.aoColumns[ iCol ];
		var oClasses = oSettings.oClasses;
		var th = $(oCol.nTh);
	
		// Try to get width information from the DOM. We can't get it from CSS
		// as we'd need to parse the CSS stylesheet. `width` option can override
		if ( ! oCol.sWidthOrig ) {
			// Width attribute
			oCol.sWidthOrig = th.attr('width') || null;
	
			// Style attribute
			var t = (th.attr('style') || '').match(/width:\s*(\d+[pxem%]+)/);
			if ( t ) {
				oCol.sWidthOrig = t[1];
			}
		}
	
		/* User specified column options */
		if ( oOptions !== undefined && oOptions !== null )
		{
			// Backwards compatibility
			_fnCompatCols( oOptions );
	
			// Map camel case parameters to their Hungarian counterparts
			_fnCamelToHungarian( DataTable.defaults.column, oOptions, true );
	
			/* Backwards compatibility for mDataProp */
			if ( oOptions.mDataProp !== undefined && !oOptions.mData )
			{
				oOptions.mData = oOptions.mDataProp;
			}
	
			if ( oOptions.sType )
			{
				oCol._sManualType = oOptions.sType;
			}
	
			// `class` is a reserved word in Javascript, so we need to provide
			// the ability to use a valid name for the camel case input
			if ( oOptions.className && ! oOptions.sClass )
			{
				oOptions.sClass = oOptions.className;
			}
			if ( oOptions.sClass ) {
				th.addClass( oOptions.sClass );
			}
	
			var origClass = oCol.sClass;
	
			$.extend( oCol, oOptions );
			_fnMap( oCol, oOptions, "sWidth", "sWidthOrig" );
	
			// Merge class from previously defined classes with this one, rather than just
			// overwriting it in the extend above
			if (origClass !== oCol.sClass) {
				oCol.sClass = origClass + ' ' + oCol.sClass;
			}
	
			/* iDataSort to be applied (backwards compatibility), but aDataSort will take
			 * priority if defined
			 */
			if ( oOptions.iDataSort !== undefined )
			{
				oCol.aDataSort = [ oOptions.iDataSort ];
			}
			_fnMap( oCol, oOptions, "aDataSort" );
		}
	
		/* Cache the data get and set functions for speed */
		var mDataSrc = oCol.mData;
		var mData = _fnGetObjectDataFn( mDataSrc );
		var mRender = oCol.mRender ? _fnGetObjectDataFn( oCol.mRender ) : null;
	
		var attrTest = function( src ) {
			return typeof src === 'string' && src.indexOf('@') !== -1;
		};
		oCol._bAttrSrc = $.isPlainObject( mDataSrc ) && (
			attrTest(mDataSrc.sort) || attrTest(mDataSrc.type) || attrTest(mDataSrc.filter)
		);
		oCol._setter = null;
	
		oCol.fnGetData = function (rowData, type, meta) {
			var innerData = mData( rowData, type, undefined, meta );
	
			return mRender && type ?
				mRender( innerData, type, rowData, meta ) :
				innerData;
		};
		oCol.fnSetData = function ( rowData, val, meta ) {
			return _fnSetObjectDataFn( mDataSrc )( rowData, val, meta );
		};
	
		// Indicate if DataTables should read DOM data as an object or array
		// Used in _fnGetRowElements
		if ( typeof mDataSrc !== 'number' && ! oCol._isArrayHost ) {
			oSettings._rowReadObject = true;
		}
	
		/* Feature sorting overrides column specific when off */
		if ( !oSettings.oFeatures.bSort )
		{
			oCol.bSortable = false;
			th.addClass( oClasses.sSortableNone ); // Have to add class here as order event isn't called
		}
	
		/* Check that the class assignment is correct for sorting */
		var bAsc = $.inArray('asc', oCol.asSorting) !== -1;
		var bDesc = $.inArray('desc', oCol.asSorting) !== -1;
		if ( !oCol.bSortable || (!bAsc && !bDesc) )
		{
			oCol.sSortingClass = oClasses.sSortableNone;
			oCol.sSortingClassJUI = "";
		}
		else if ( bAsc && !bDesc )
		{
			oCol.sSortingClass = oClasses.sSortableAsc;
			oCol.sSortingClassJUI = oClasses.sSortJUIAscAllowed;
		}
		else if ( !bAsc && bDesc )
		{
			oCol.sSortingClass = oClasses.sSortableDesc;
			oCol.sSortingClassJUI = oClasses.sSortJUIDescAllowed;
		}
		else
		{
			oCol.sSortingClass = oClasses.sSortable;
			oCol.sSortingClassJUI = oClasses.sSortJUI;
		}
	}
	
	
	/**
	 * Adjust the table column widths for new data. Note: you would probably want to
	 * do a redraw after calling this function!
	 *  @param {object} settings dataTables settings object
	 *  @memberof DataTable#oApi
	 */
	function _fnAdjustColumnSizing ( settings )
	{
		/* Not interested in doing column width calculation if auto-width is disabled */
		if ( settings.oFeatures.bAutoWidth !== false )
		{
			var columns = settings.aoColumns;
	
			_fnCalculateColumnWidths( settings );
			for ( var i=0 , iLen=columns.length ; i<iLen ; i++ )
			{
				columns[i].nTh.style.width = columns[i].sWidth;
			}
		}
	
		var scroll = settings.oScroll;
		if ( scroll.sY !== '' || scroll.sX !== '')
		{
			_fnScrollDraw( settings );
		}
	
		_fnCallbackFire( settings, null, 'column-sizing', [settings] );
	}
	
	
	/**
	 * Convert the index of a visible column to the index in the data array (take account
	 * of hidden columns)
	 *  @param {object} oSettings dataTables settings object
	 *  @param {int} iMatch Visible column index to lookup
	 *  @returns {int} i the data index
	 *  @memberof DataTable#oApi
	 */
	function _fnVisibleToColumnIndex( oSettings, iMatch )
	{
		var aiVis = _fnGetColumns( oSettings, 'bVisible' );
	
		return typeof aiVis[iMatch] === 'number' ?
			aiVis[iMatch] :
			null;
	}
	
	
	/**
	 * Convert the index of an index in the data array and convert it to the visible
	 *   column index (take account of hidden columns)
	 *  @param {int} iMatch Column index to lookup
	 *  @param {object} oSettings dataTables settings object
	 *  @returns {int} i the data index
	 *  @memberof DataTable#oApi
	 */
	function _fnColumnIndexToVisible( oSettings, iMatch )
	{
		var aiVis = _fnGetColumns( oSettings, 'bVisible' );
		var iPos = $.inArray( iMatch, aiVis );
	
		return iPos !== -1 ? iPos : null;
	}
	
	
	/**
	 * Get the number of visible columns
	 *  @param {object} oSettings dataTables settings object
	 *  @returns {int} i the number of visible columns
	 *  @memberof DataTable#oApi
	 */
	function _fnVisbleColumns( oSettings )
	{
		var vis = 0;
	
		// No reduce in IE8, use a loop for now
		$.each( oSettings.aoColumns, function ( i, col ) {
			if ( col.bVisible && $(col.nTh).css('display') !== 'none' ) {
				vis++;
			}
		} );
	
		return vis;
	}
	
	
	/**
	 * Get an array of column indexes that match a given property
	 *  @param {object} oSettings dataTables settings object
	 *  @param {string} sParam Parameter in aoColumns to look for - typically
	 *    bVisible or bSearchable
	 *  @returns {array} Array of indexes with matched properties
	 *  @memberof DataTable#oApi
	 */
	function _fnGetColumns( oSettings, sParam )
	{
		var a = [];
	
		$.map( oSettings.aoColumns, function(val, i) {
			if ( val[sParam] ) {
				a.push( i );
			}
		} );
	
		return a;
	}
	
	
	/**
	 * Calculate the 'type' of a column
	 *  @param {object} settings dataTables settings object
	 *  @memberof DataTable#oApi
	 */
	function _fnColumnTypes ( settings )
	{
		var columns = settings.aoColumns;
		var data = settings.aoData;
		var types = DataTable.ext.type.detect;
		var i, ien, j, jen, k, ken;
		var col, cell, detectedType, cache;
	
		// For each column, spin over the 
		for ( i=0, ien=columns.length ; i<ien ; i++ ) {
			col = columns[i];
			cache = [];
	
			if ( ! col.sType && col._sManualType ) {
				col.sType = col._sManualType;
			}
			else if ( ! col.sType ) {
				for ( j=0, jen=types.length ; j<jen ; j++ ) {
					for ( k=0, ken=data.length ; k<ken ; k++ ) {
						// Use a cache array so we only need to get the type data
						// from the formatter once (when using multiple detectors)
						if ( cache[k] === undefined ) {
							cache[k] = _fnGetCellData( settings, k, i, 'type' );
						}
	
						detectedType = types[j]( cache[k], settings );
	
						// If null, then this type can't apply to this column, so
						// rather than testing all cells, break out. There is an
						// exception for the last type which is `html`. We need to
						// scan all rows since it is possible to mix string and HTML
						// types
						if ( ! detectedType && j !== types.length-1 ) {
							break;
						}
	
						// Only a single match is needed for html type since it is
						// bottom of the pile and very similar to string - but it
						// must not be empty
						if ( detectedType === 'html' && ! _empty(cache[k]) ) {
							break;
						}
					}
	
					// Type is valid for all data points in the column - use this
					// type
					if ( detectedType ) {
						col.sType = detectedType;
						break;
					}
				}
	
				// Fall back - if no type was detected, always use string
				if ( ! col.sType ) {
					col.sType = 'string';
				}
			}
		}
	}
	
	
	/**
	 * Take the column definitions and static columns arrays and calculate how
	 * they relate to column indexes. The callback function will then apply the
	 * definition found for a column to a suitable configuration object.
	 *  @param {object} oSettings dataTables settings object
	 *  @param {array} aoColDefs The aoColumnDefs array that is to be applied
	 *  @param {array} aoCols The aoColumns array that defines columns individually
	 *  @param {function} fn Callback function - takes two parameters, the calculated
	 *    column index and the definition for that column.
	 *  @memberof DataTable#oApi
	 */
	function _fnApplyColumnDefs( oSettings, aoColDefs, aoCols, fn )
	{
		var i, iLen, j, jLen, k, kLen, def;
		var columns = oSettings.aoColumns;
	
		// Column definitions with aTargets
		if ( aoColDefs )
		{
			/* Loop over the definitions array - loop in reverse so first instance has priority */
			for ( i=aoColDefs.length-1 ; i>=0 ; i-- )
			{
				def = aoColDefs[i];
	
				/* Each definition can target multiple columns, as it is an array */
				var aTargets = def.target !== undefined
					? def.target
					: def.targets !== undefined
						? def.targets
						: def.aTargets;
	
				if ( ! Array.isArray( aTargets ) )
				{
					aTargets = [ aTargets ];
				}
	
				for ( j=0, jLen=aTargets.length ; j<jLen ; j++ )
				{
					if ( typeof aTargets[j] === 'number' && aTargets[j] >= 0 )
					{
						/* Add columns that we don't yet know about */
						while( columns.length <= aTargets[j] )
						{
							_fnAddColumn( oSettings );
						}
	
						/* Integer, basic index */
						fn( aTargets[j], def );
					}
					else if ( typeof aTargets[j] === 'number' && aTargets[j] < 0 )
					{
						/* Negative integer, right to left column counting */
						fn( columns.length+aTargets[j], def );
					}
					else if ( typeof aTargets[j] === 'string' )
					{
						/* Class name matching on TH element */
						for ( k=0, kLen=columns.length ; k<kLen ; k++ )
						{
							if ( aTargets[j] == "_all" ||
							     $(columns[k].nTh).hasClass( aTargets[j] ) )
							{
								fn( k, def );
							}
						}
					}
				}
			}
		}
	
		// Statically defined columns array
		if ( aoCols )
		{
			for ( i=0, iLen=aoCols.length ; i<iLen ; i++ )
			{
				fn( i, aoCols[i] );
			}
		}
	}
	
	/**
	 * Add a data array to the table, creating DOM node etc. This is the parallel to
	 * _fnGatherData, but for adding rows from a Javascript source, rather than a
	 * DOM source.
	 *  @param {object} oSettings dataTables settings object
	 *  @param {array} aData data array to be added
	 *  @param {node} [nTr] TR element to add to the table - optional. If not given,
	 *    DataTables will create a row automatically
	 *  @param {array} [anTds] Array of TD|TH elements for the row - must be given
	 *    if nTr is.
	 *  @returns {int} >=0 if successful (index of new aoData entry), -1 if failed
	 *  @memberof DataTable#oApi
	 */
	function _fnAddData ( oSettings, aDataIn, nTr, anTds )
	{
		/* Create the object for storing information about this new row */
		var iRow = oSettings.aoData.length;
		var oData = $.extend( true, {}, DataTable.models.oRow, {
			src: nTr ? 'dom' : 'data',
			idx: iRow
		} );
	
		oData._aData = aDataIn;
		oSettings.aoData.push( oData );
	
		/* Create the cells */
		var nTd, sThisType;
		var columns = oSettings.aoColumns;
	
		// Invalidate the column types as the new data needs to be revalidated
		for ( var i=0, iLen=columns.length ; i<iLen ; i++ )
		{
			columns[i].sType = null;
		}
	
		/* Add to the display array */
		oSettings.aiDisplayMaster.push( iRow );
	
		var id = oSettings.rowIdFn( aDataIn );
		if ( id !== undefined ) {
			oSettings.aIds[ id ] = oData;
		}
	
		/* Create the DOM information, or register it if already present */
		if ( nTr || ! oSettings.oFeatures.bDeferRender )
		{
			_fnCreateTr( oSettings, iRow, nTr, anTds );
		}
	
		return iRow;
	}
	
	
	/**
	 * Add one or more TR elements to the table. Generally we'd expect to
	 * use this for reading data from a DOM sourced table, but it could be
	 * used for an TR element. Note that if a TR is given, it is used (i.e.
	 * it is not cloned).
	 *  @param {object} settings dataTables settings object
	 *  @param {array|node|jQuery} trs The TR element(s) to add to the table
	 *  @returns {array} Array of indexes for the added rows
	 *  @memberof DataTable#oApi
	 */
	function _fnAddTr( settings, trs )
	{
		var row;
	
		// Allow an individual node to be passed in
		if ( ! (trs instanceof $) ) {
			trs = $(trs);
		}
	
		return trs.map( function (i, el) {
			row = _fnGetRowElements( settings, el );
			return _fnAddData( settings, row.data, el, row.cells );
		} );
	}
	
	
	/**
	 * Take a TR element and convert it to an index in aoData
	 *  @param {object} oSettings dataTables settings object
	 *  @param {node} n the TR element to find
	 *  @returns {int} index if the node is found, null if not
	 *  @memberof DataTable#oApi
	 */
	function _fnNodeToDataIndex( oSettings, n )
	{
		return (n._DT_RowIndex!==undefined) ? n._DT_RowIndex : null;
	}
	
	
	/**
	 * Take a TD element and convert it into a column data index (not the visible index)
	 *  @param {object} oSettings dataTables settings object
	 *  @param {int} iRow The row number the TD/TH can be found in
	 *  @param {node} n The TD/TH element to find
	 *  @returns {int} index if the node is found, -1 if not
	 *  @memberof DataTable#oApi
	 */
	function _fnNodeToColumnIndex( oSettings, iRow, n )
	{
		return $.inArray( n, oSettings.aoData[ iRow ].anCells );
	}
	
	
	/**
	 * Get the data for a given cell from the internal cache, taking into account data mapping
	 *  @param {object} settings dataTables settings object
	 *  @param {int} rowIdx aoData row id
	 *  @param {int} colIdx Column index
	 *  @param {string} type data get type ('display', 'type' 'filter|search' 'sort|order')
	 *  @returns {*} Cell data
	 *  @memberof DataTable#oApi
	 */
	function _fnGetCellData( settings, rowIdx, colIdx, type )
	{
		if (type === 'search') {
			type = 'filter';
		}
		else if (type === 'order') {
			type = 'sort';
		}
	
		var draw           = settings.iDraw;
		var col            = settings.aoColumns[colIdx];
		var rowData        = settings.aoData[rowIdx]._aData;
		var defaultContent = col.sDefaultContent;
		var cellData       = col.fnGetData( rowData, type, {
			settings: settings,
			row:      rowIdx,
			col:      colIdx
		} );
	
		if ( cellData === undefined ) {
			if ( settings.iDrawError != draw && defaultContent === null ) {
				_fnLog( settings, 0, "Requested unknown parameter "+
					(typeof col.mData=='function' ? '{function}' : "'"+col.mData+"'")+
					" for row "+rowIdx+", column "+colIdx, 4 );
				settings.iDrawError = draw;
			}
			return defaultContent;
		}
	
		// When the data source is null and a specific data type is requested (i.e.
		// not the original data), we can use default column data
		if ( (cellData === rowData || cellData === null) && defaultContent !== null && type !== undefined ) {
			cellData = defaultContent;
		}
		else if ( typeof cellData === 'function' ) {
			// If the data source is a function, then we run it and use the return,
			// executing in the scope of the data object (for instances)
			return cellData.call( rowData );
		}
	
		if ( cellData === null && type === 'display' ) {
			return '';
		}
	
		if ( type === 'filter' ) {
			var fomatters = DataTable.ext.type.search;
	
			if ( fomatters[ col.sType ] ) {
				cellData = fomatters[ col.sType ]( cellData );
			}
		}
	
		return cellData;
	}
	
	
	/**
	 * Set the value for a specific cell, into the internal data cache
	 *  @param {object} settings dataTables settings object
	 *  @param {int} rowIdx aoData row id
	 *  @param {int} colIdx Column index
	 *  @param {*} val Value to set
	 *  @memberof DataTable#oApi
	 */
	function _fnSetCellData( settings, rowIdx, colIdx, val )
	{
		var col     = settings.aoColumns[colIdx];
		var rowData = settings.aoData[rowIdx]._aData;
	
		col.fnSetData( rowData, val, {
			settings: settings,
			row:      rowIdx,
			col:      colIdx
		}  );
	}
	
	
	// Private variable that is used to match action syntax in the data property object
	var __reArray = /\[.*?\]$/;
	var __reFn = /\(\)$/;
	
	/**
	 * Split string on periods, taking into account escaped periods
	 * @param  {string} str String to split
	 * @return {array} Split string
	 */
	function _fnSplitObjNotation( str )
	{
		return $.map( str.match(/(\\.|[^\.])+/g) || [''], function ( s ) {
			return s.replace(/\\\./g, '.');
		} );
	}
	
	
	/**
	 * Return a function that can be used to get data from a source object, taking
	 * into account the ability to use nested objects as a source
	 *  @param {string|int|function} mSource The data source for the object
	 *  @returns {function} Data get function
	 *  @memberof DataTable#oApi
	 */
	var _fnGetObjectDataFn = DataTable.util.get;
	
	
	/**
	 * Return a function that can be used to set data from a source object, taking
	 * into account the ability to use nested objects as a source
	 *  @param {string|int|function} mSource The data source for the object
	 *  @returns {function} Data set function
	 *  @memberof DataTable#oApi
	 */
	var _fnSetObjectDataFn = DataTable.util.set;
	
	
	/**
	 * Return an array with the full table data
	 *  @param {object} oSettings dataTables settings object
	 *  @returns array {array} aData Master data array
	 *  @memberof DataTable#oApi
	 */
	function _fnGetDataMaster ( settings )
	{
		return _pluck( settings.aoData, '_aData' );
	}
	
	
	/**
	 * Nuke the table
	 *  @param {object} oSettings dataTables settings object
	 *  @memberof DataTable#oApi
	 */
	function _fnClearTable( settings )
	{
		settings.aoData.length = 0;
		settings.aiDisplayMaster.length = 0;
		settings.aiDisplay.length = 0;
		settings.aIds = {};
	}
	
	
	 /**
	 * Take an array of integers (index array) and remove a target integer (value - not
	 * the key!)
	 *  @param {array} a Index array to target
	 *  @param {int} iTarget value to find
	 *  @memberof DataTable#oApi
	 */
	function _fnDeleteIndex( a, iTarget, splice )
	{
		var iTargetIndex = -1;
	
		for ( var i=0, iLen=a.length ; i<iLen ; i++ )
		{
			if ( a[i] == iTarget )
			{
				iTargetIndex = i;
			}
			else if ( a[i] > iTarget )
			{
				a[i]--;
			}
		}
	
		if ( iTargetIndex != -1 && splice === undefined )
		{
			a.splice( iTargetIndex, 1 );
		}
	}
	
	
	/**
	 * Mark cached data as invalid such that a re-read of the data will occur when
	 * the cached data is next requested. Also update from the data source object.
	 *
	 * @param {object} settings DataTables settings object
	 * @param {int}    rowIdx   Row index to invalidate
	 * @param {string} [src]    Source to invalidate from: undefined, 'auto', 'dom'
	 *     or 'data'
	 * @param {int}    [colIdx] Column index to invalidate. If undefined the whole
	 *     row will be invalidated
	 * @memberof DataTable#oApi
	 *
	 * @todo For the modularisation of v1.11 this will need to become a callback, so
	 *   the sort and filter methods can subscribe to it. That will required
	 *   initialisation options for sorting, which is why it is not already baked in
	 */
	function _fnInvalidate( settings, rowIdx, src, colIdx )
	{
		var row = settings.aoData[ rowIdx ];
		var i, ien;
		var cellWrite = function ( cell, col ) {
			// This is very frustrating, but in IE if you just write directly
			// to innerHTML, and elements that are overwritten are GC'ed,
			// even if there is a reference to them elsewhere
			while ( cell.childNodes.length ) {
				cell.removeChild( cell.firstChild );
			}
	
			cell.innerHTML = _fnGetCellData( settings, rowIdx, col, 'display' );
		};
	
		// Are we reading last data from DOM or the data object?
		if ( src === 'dom' || ((! src || src === 'auto') && row.src === 'dom') ) {
			// Read the data from the DOM
			row._aData = _fnGetRowElements(
					settings, row, colIdx, colIdx === undefined ? undefined : row._aData
				)
				.data;
		}
		else {
			// Reading from data object, update the DOM
			var cells = row.anCells;
	
			if ( cells ) {
				if ( colIdx !== undefined ) {
					cellWrite( cells[colIdx], colIdx );
				}
				else {
					for ( i=0, ien=cells.length ; i<ien ; i++ ) {
						cellWrite( cells[i], i );
					}
				}
			}
		}
	
		// For both row and cell invalidation, the cached data for sorting and
		// filtering is nulled out
		row._aSortData = null;
		row._aFilterData = null;
	
		// Invalidate the type for a specific column (if given) or all columns since
		// the data might have changed
		var cols = settings.aoColumns;
		if ( colIdx !== undefined ) {
			cols[ colIdx ].sType = null;
		}
		else {
			for ( i=0, ien=cols.length ; i<ien ; i++ ) {
				cols[i].sType = null;
			}
	
			// Update DataTables special `DT_*` attributes for the row
			_fnRowAttributes( settings, row );
		}
	}
	
	
	/**
	 * Build a data source object from an HTML row, reading the contents of the
	 * cells that are in the row.
	 *
	 * @param {object} settings DataTables settings object
	 * @param {node|object} TR element from which to read data or existing row
	 *   object from which to re-read the data from the cells
	 * @param {int} [colIdx] Optional column index
	 * @param {array|object} [d] Data source object. If `colIdx` is given then this
	 *   parameter should also be given and will be used to write the data into.
	 *   Only the column in question will be written
	 * @returns {object} Object with two parameters: `data` the data read, in
	 *   document order, and `cells` and array of nodes (they can be useful to the
	 *   caller, so rather than needing a second traversal to get them, just return
	 *   them from here).
	 * @memberof DataTable#oApi
	 */
	function _fnGetRowElements( settings, row, colIdx, d )
	{
		var
			tds = [],
			td = row.firstChild,
			name, col, o, i=0, contents,
			columns = settings.aoColumns,
			objectRead = settings._rowReadObject;
	
		// Allow the data object to be passed in, or construct
		d = d !== undefined ?
			d :
			objectRead ?
				{} :
				[];
	
		var attr = function ( str, td  ) {
			if ( typeof str === 'string' ) {
				var idx = str.indexOf('@');
	
				if ( idx !== -1 ) {
					var attr = str.substring( idx+1 );
					var setter = _fnSetObjectDataFn( str );
					setter( d, td.getAttribute( attr ) );
				}
			}
		};
	
		// Read data from a cell and store into the data object
		var cellProcess = function ( cell ) {
			if ( colIdx === undefined || colIdx === i ) {
				col = columns[i];
				contents = (cell.innerHTML).trim();
	
				if ( col && col._bAttrSrc ) {
					var setter = _fnSetObjectDataFn( col.mData._ );
					setter( d, contents );
	
					attr( col.mData.sort, cell );
					attr( col.mData.type, cell );
					attr( col.mData.filter, cell );
				}
				else {
					// Depending on the `data` option for the columns the data can
					// be read to either an object or an array.
					if ( objectRead ) {
						if ( ! col._setter ) {
							// Cache the setter function
							col._setter = _fnSetObjectDataFn( col.mData );
						}
						col._setter( d, contents );
					}
					else {
						d[i] = contents;
					}
				}
			}
	
			i++;
		};
	
		if ( td ) {
			// `tr` element was passed in
			while ( td ) {
				name = td.nodeName.toUpperCase();
	
				if ( name == "TD" || name == "TH" ) {
					cellProcess( td );
					tds.push( td );
				}
	
				td = td.nextSibling;
			}
		}
		else {
			// Existing row object passed in
			tds = row.anCells;
	
			for ( var j=0, jen=tds.length ; j<jen ; j++ ) {
				cellProcess( tds[j] );
			}
		}
	
		// Read the ID from the DOM if present
		var rowNode = row.firstChild ? row : row.nTr;
	
		if ( rowNode ) {
			var id = rowNode.getAttribute( 'id' );
	
			if ( id ) {
				_fnSetObjectDataFn( settings.rowId )( d, id );
			}
		}
	
		return {
			data: d,
			cells: tds
		};
	}
	/**
	 * Create a new TR element (and it's TD children) for a row
	 *  @param {object} oSettings dataTables settings object
	 *  @param {int} iRow Row to consider
	 *  @param {node} [nTrIn] TR element to add to the table - optional. If not given,
	 *    DataTables will create a row automatically
	 *  @param {array} [anTds] Array of TD|TH elements for the row - must be given
	 *    if nTr is.
	 *  @memberof DataTable#oApi
	 */
	function _fnCreateTr ( oSettings, iRow, nTrIn, anTds )
	{
		var
			row = oSettings.aoData[iRow],
			rowData = row._aData,
			cells = [],
			nTr, nTd, oCol,
			i, iLen, create;
	
		if ( row.nTr === null )
		{
			nTr = nTrIn || document.createElement('tr');
	
			row.nTr = nTr;
			row.anCells = cells;
	
			/* Use a private property on the node to allow reserve mapping from the node
			 * to the aoData array for fast look up
			 */
			nTr._DT_RowIndex = iRow;
	
			/* Special parameters can be given by the data source to be used on the row */
			_fnRowAttributes( oSettings, row );
	
			/* Process each column */
			for ( i=0, iLen=oSettings.aoColumns.length ; i<iLen ; i++ )
			{
				oCol = oSettings.aoColumns[i];
				create = nTrIn ? false : true;
	
				nTd = create ? document.createElement( oCol.sCellType ) : anTds[i];
	
				if (! nTd) {
					_fnLog( oSettings, 0, 'Incorrect column count', 18 );
				}
	
				nTd._DT_CellIndex = {
					row: iRow,
					column: i
				};
				
				cells.push( nTd );
	
				// Need to create the HTML if new, or if a rendering function is defined
				if ( create || ((oCol.mRender || oCol.mData !== i) &&
					 (!$.isPlainObject(oCol.mData) || oCol.mData._ !== i+'.display')
				)) {
					nTd.innerHTML = _fnGetCellData( oSettings, iRow, i, 'display' );
				}
	
				/* Add user defined class */
				if ( oCol.sClass )
				{
					nTd.className += ' '+oCol.sClass;
				}
	
				// Visibility - add or remove as required
				if ( oCol.bVisible && ! nTrIn )
				{
					nTr.appendChild( nTd );
				}
				else if ( ! oCol.bVisible && nTrIn )
				{
					nTd.parentNode.removeChild( nTd );
				}
	
				if ( oCol.fnCreatedCell )
				{
					oCol.fnCreatedCell.call( oSettings.oInstance,
						nTd, _fnGetCellData( oSettings, iRow, i ), rowData, iRow, i
					);
				}
			}
	
			_fnCallbackFire( oSettings, 'aoRowCreatedCallback', null, [nTr, rowData, iRow, cells] );
		}
	}
	
	
	/**
	 * Add attributes to a row based on the special `DT_*` parameters in a data
	 * source object.
	 *  @param {object} settings DataTables settings object
	 *  @param {object} DataTables row object for the row to be modified
	 *  @memberof DataTable#oApi
	 */
	function _fnRowAttributes( settings, row )
	{
		var tr = row.nTr;
		var data = row._aData;
	
		if ( tr ) {
			var id = settings.rowIdFn( data );
	
			if ( id ) {
				tr.id = id;
			}
	
			if ( data.DT_RowClass ) {
				// Remove any classes added by DT_RowClass before
				var a = data.DT_RowClass.split(' ');
				row.__rowc = row.__rowc ?
					_unique( row.__rowc.concat( a ) ) :
					a;
	
				$(tr)
					.removeClass( row.__rowc.join(' ') )
					.addClass( data.DT_RowClass );
			}
	
			if ( data.DT_RowAttr ) {
				$(tr).attr( data.DT_RowAttr );
			}
	
			if ( data.DT_RowData ) {
				$(tr).data( data.DT_RowData );
			}
		}
	}
	
	
	/**
	 * Create the HTML header for the table
	 *  @param {object} oSettings dataTables settings object
	 *  @memberof DataTable#oApi
	 */
	function _fnBuildHead( oSettings )
	{
		var i, ien, cell, row, column;
		var thead = oSettings.nTHead;
		var tfoot = oSettings.nTFoot;
		var createHeader = $('th, td', thead).length === 0;
		var classes = oSettings.oClasses;
		var columns = oSettings.aoColumns;
	
		if ( createHeader ) {
			row = $('<tr/>').appendTo( thead );
		}
	
		for ( i=0, ien=columns.length ; i<ien ; i++ ) {
			column = columns[i];
			cell = $( column.nTh ).addClass( column.sClass );
	
			if ( createHeader ) {
				cell.appendTo( row );
			}
	
			// 1.11 move into sorting
			if ( oSettings.oFeatures.bSort ) {
				cell.addClass( column.sSortingClass );
	
				if ( column.bSortable !== false ) {
					cell
						.attr( 'tabindex', oSettings.iTabIndex )
						.attr( 'aria-controls', oSettings.sTableId );
	
					_fnSortAttachListener( oSettings, column.nTh, i );
				}
			}
	
			if ( column.sTitle != cell[0].innerHTML ) {
				cell.html( column.sTitle );
			}
	
			_fnRenderer( oSettings, 'header' )(
				oSettings, cell, column, classes
			);
		}
	
		if ( createHeader ) {
			_fnDetectHeader( oSettings.aoHeader, thead );
		}
	
		/* Deal with the footer - add classes if required */
		$(thead).children('tr').children('th, td').addClass( classes.sHeaderTH );
		$(tfoot).children('tr').children('th, td').addClass( classes.sFooterTH );
	
		// Cache the footer cells. Note that we only take the cells from the first
		// row in the footer. If there is more than one row the user wants to
		// interact with, they need to use the table().foot() method. Note also this
		// allows cells to be used for multiple columns using colspan
		if ( tfoot !== null ) {
			var cells = oSettings.aoFooter[0];
	
			for ( i=0, ien=cells.length ; i<ien ; i++ ) {
				column = columns[i];
	
				if (column) {
					column.nTf = cells[i].cell;
		
					if ( column.sClass ) {
						$(column.nTf).addClass( column.sClass );
					}
				}
				else {
					_fnLog( oSettings, 0, 'Incorrect column count', 18 );
				}
			}
		}
	}
	
	
	/**
	 * Draw the header (or footer) element based on the column visibility states. The
	 * methodology here is to use the layout array from _fnDetectHeader, modified for
	 * the instantaneous column visibility, to construct the new layout. The grid is
	 * traversed over cell at a time in a rows x columns grid fashion, although each
	 * cell insert can cover multiple elements in the grid - which is tracks using the
	 * aApplied array. Cell inserts in the grid will only occur where there isn't
	 * already a cell in that position.
	 *  @param {object} oSettings dataTables settings object
	 *  @param array {objects} aoSource Layout array from _fnDetectHeader
	 *  @param {boolean} [bIncludeHidden=false] If true then include the hidden columns in the calc,
	 *  @memberof DataTable#oApi
	 */
	function _fnDrawHead( oSettings, aoSource, bIncludeHidden )
	{
		var i, iLen, j, jLen, k, kLen, n, nLocalTr;
		var aoLocal = [];
		var aApplied = [];
		var iColumns = oSettings.aoColumns.length;
		var iRowspan, iColspan;
	
		if ( ! aoSource )
		{
			return;
		}
	
		if (  bIncludeHidden === undefined )
		{
			bIncludeHidden = false;
		}
	
		/* Make a copy of the master layout array, but without the visible columns in it */
		for ( i=0, iLen=aoSource.length ; i<iLen ; i++ )
		{
			aoLocal[i] = aoSource[i].slice();
			aoLocal[i].nTr = aoSource[i].nTr;
	
			/* Remove any columns which are currently hidden */
			for ( j=iColumns-1 ; j>=0 ; j-- )
			{
				if ( !oSettings.aoColumns[j].bVisible && !bIncludeHidden )
				{
					aoLocal[i].splice( j, 1 );
				}
			}
	
			/* Prep the applied array - it needs an element for each row */
			aApplied.push( [] );
		}
	
		for ( i=0, iLen=aoLocal.length ; i<iLen ; i++ )
		{
			nLocalTr = aoLocal[i].nTr;
	
			/* All cells are going to be replaced, so empty out the row */
			if ( nLocalTr )
			{
				while( (n = nLocalTr.firstChild) )
				{
					nLocalTr.removeChild( n );
				}
			}
	
			for ( j=0, jLen=aoLocal[i].length ; j<jLen ; j++ )
			{
				iRowspan = 1;
				iColspan = 1;
	
				/* Check to see if there is already a cell (row/colspan) covering our target
				 * insert point. If there is, then there is nothing to do.
				 */
				if ( aApplied[i][j] === undefined )
				{
					nLocalTr.appendChild( aoLocal[i][j].cell );
					aApplied[i][j] = 1;
	
					/* Expand the cell to cover as many rows as needed */
					while ( aoLocal[i+iRowspan] !== undefined &&
					        aoLocal[i][j].cell == aoLocal[i+iRowspan][j].cell )
					{
						aApplied[i+iRowspan][j] = 1;
						iRowspan++;
					}
	
					/* Expand the cell to cover as many columns as needed */
					while ( aoLocal[i][j+iColspan] !== undefined &&
					        aoLocal[i][j].cell == aoLocal[i][j+iColspan].cell )
					{
						/* Must update the applied array over the rows for the columns */
						for ( k=0 ; k<iRowspan ; k++ )
						{
							aApplied[i+k][j+iColspan] = 1;
						}
						iColspan++;
					}
	
					/* Do the actual expansion in the DOM */
					$(aoLocal[i][j].cell)
						.attr('rowspan', iRowspan)
						.attr('colspan', iColspan);
				}
			}
		}
	}
	
	
	/**
	 * Insert the required TR nodes into the table for display
	 *  @param {object} oSettings dataTables settings object
	 *  @param ajaxComplete true after ajax call to complete rendering
	 *  @memberof DataTable#oApi
	 */
	function _fnDraw( oSettings, ajaxComplete )
	{
		// Allow for state saving and a custom start position
		_fnStart( oSettings );
	
		/* Provide a pre-callback function which can be used to cancel the draw is false is returned */
		var aPreDraw = _fnCallbackFire( oSettings, 'aoPreDrawCallback', 'preDraw', [oSettings] );
		if ( $.inArray( false, aPreDraw ) !== -1 )
		{
			_fnProcessingDisplay( oSettings, false );
			return;
		}
	
		var anRows = [];
		var iRowCount = 0;
		var asStripeClasses = oSettings.asStripeClasses;
		var iStripes = asStripeClasses.length;
		var oLang = oSettings.oLanguage;
		var bServerSide = _fnDataSource( oSettings ) == 'ssp';
		var aiDisplay = oSettings.aiDisplay;
		var iDisplayStart = oSettings._iDisplayStart;
		var iDisplayEnd = oSettings.fnDisplayEnd();
	
		oSettings.bDrawing = true;
	
		/* Server-side processing draw intercept */
		if ( oSettings.bDeferLoading )
		{
			oSettings.bDeferLoading = false;
			oSettings.iDraw++;
			_fnProcessingDisplay( oSettings, false );
		}
		else if ( !bServerSide )
		{
			oSettings.iDraw++;
		}
		else if ( !oSettings.bDestroying && !ajaxComplete)
		{
			_fnAjaxUpdate( oSettings );
			return;
		}
	
		if ( aiDisplay.length !== 0 )
		{
			var iStart = bServerSide ? 0 : iDisplayStart;
			var iEnd = bServerSide ? oSettings.aoData.length : iDisplayEnd;
	
			for ( var j=iStart ; j<iEnd ; j++ )
			{
				var iDataIndex = aiDisplay[j];
				var aoData = oSettings.aoData[ iDataIndex ];
				if ( aoData.nTr === null )
				{
					_fnCreateTr( oSettings, iDataIndex );
				}
	
				var nRow = aoData.nTr;
	
				/* Remove the old striping classes and then add the new one */
				if ( iStripes !== 0 )
				{
					var sStripe = asStripeClasses[ iRowCount % iStripes ];
					if ( aoData._sRowStripe != sStripe )
					{
						$(nRow).removeClass( aoData._sRowStripe ).addClass( sStripe );
						aoData._sRowStripe = sStripe;
					}
				}
	
				// Row callback functions - might want to manipulate the row
				// iRowCount and j are not currently documented. Are they at all
				// useful?
				_fnCallbackFire( oSettings, 'aoRowCallback', null,
					[nRow, aoData._aData, iRowCount, j, iDataIndex] );
	
				anRows.push( nRow );
				iRowCount++;
			}
		}
		else
		{
			/* Table is empty - create a row with an empty message in it */
			var sZero = oLang.sZeroRecords;
			if ( oSettings.iDraw == 1 &&  _fnDataSource( oSettings ) == 'ajax' )
			{
				sZero = oLang.sLoadingRecords;
			}
			else if ( oLang.sEmptyTable && oSettings.fnRecordsTotal() === 0 )
			{
				sZero = oLang.sEmptyTable;
			}
	
			anRows[ 0 ] = $( '<tr/>', { 'class': iStripes ? asStripeClasses[0] : '' } )
				.append( $('<td />', {
					'valign':  'top',
					'colSpan': _fnVisbleColumns( oSettings ),
					'class':   oSettings.oClasses.sRowEmpty
				} ).html( sZero ) )[0];
		}
	
		/* Header and footer callbacks */
		_fnCallbackFire( oSettings, 'aoHeaderCallback', 'header', [ $(oSettings.nTHead).children('tr')[0],
			_fnGetDataMaster( oSettings ), iDisplayStart, iDisplayEnd, aiDisplay ] );
	
		_fnCallbackFire( oSettings, 'aoFooterCallback', 'footer', [ $(oSettings.nTFoot).children('tr')[0],
			_fnGetDataMaster( oSettings ), iDisplayStart, iDisplayEnd, aiDisplay ] );
	
		var body = $(oSettings.nTBody);
	
		body.children().detach();
		body.append( $(anRows) );
	
		/* Call all required callback functions for the end of a draw */
		_fnCallbackFire( oSettings, 'aoDrawCallback', 'draw', [oSettings] );
	
		/* Draw is complete, sorting and filtering must be as well */
		oSettings.bSorted = false;
		oSettings.bFiltered = false;
		oSettings.bDrawing = false;
	}
	
	
	/**
	 * Redraw the table - taking account of the various features which are enabled
	 *  @param {object} oSettings dataTables settings object
	 *  @param {boolean} [holdPosition] Keep the current paging position. By default
	 *    the paging is reset to the first page
	 *  @memberof DataTable#oApi
	 */
	function _fnReDraw( settings, holdPosition )
	{
		var
			features = settings.oFeatures,
			sort     = features.bSort,
			filter   = features.bFilter;
	
		if ( sort ) {
			_fnSort( settings );
		}
	
		if ( filter ) {
			_fnFilterComplete( settings, settings.oPreviousSearch );
		}
		else {
			// No filtering, so we want to just use the display master
			settings.aiDisplay = settings.aiDisplayMaster.slice();
		}
	
		if ( holdPosition !== true ) {
			settings._iDisplayStart = 0;
		}
	
		// Let any modules know about the draw hold position state (used by
		// scrolling internally)
		settings._drawHold = holdPosition;
	
		_fnDraw( settings );
	
		settings._drawHold = false;
	}
	
	
	/**
	 * Add the options to the page HTML for the table
	 *  @param {object} oSettings dataTables settings object
	 *  @memberof DataTable#oApi
	 */
	function _fnAddOptionsHtml ( oSettings )
	{
		var classes = oSettings.oClasses;
		var table = $(oSettings.nTable);
		var holding = $('<div/>').insertBefore( table ); // Holding element for speed
		var features = oSettings.oFeatures;
	
		// All DataTables are wrapped in a div
		var insert = $('<div/>', {
			id:      oSettings.sTableId+'_wrapper',
			'class': classes.sWrapper + (oSettings.nTFoot ? '' : ' '+classes.sNoFooter)
		} );
	
		oSettings.nHolding = holding[0];
		oSettings.nTableWrapper = insert[0];
		oSettings.nTableReinsertBefore = oSettings.nTable.nextSibling;
	
		/* Loop over the user set positioning and place the elements as needed */
		var aDom = oSettings.sDom.split('');
		var featureNode, cOption, nNewNode, cNext, sAttr, j;
		for ( var i=0 ; i<aDom.length ; i++ )
		{
			featureNode = null;
			cOption = aDom[i];
	
			if ( cOption == '<' )
			{
				/* New container div */
				nNewNode = $('<div/>')[0];
	
				/* Check to see if we should append an id and/or a class name to the container */
				cNext = aDom[i+1];
				if ( cNext == "'" || cNext == '"' )
				{
					sAttr = "";
					j = 2;
					while ( aDom[i+j] != cNext )
					{
						sAttr += aDom[i+j];
						j++;
					}
	
					/* Replace jQuery UI constants @todo depreciated */
					if ( sAttr == "H" )
					{
						sAttr = classes.sJUIHeader;
					}
					else if ( sAttr == "F" )
					{
						sAttr = classes.sJUIFooter;
					}
	
					/* The attribute can be in the format of "#id.class", "#id" or "class" This logic
					 * breaks the string into parts and applies them as needed
					 */
					if ( sAttr.indexOf('.') != -1 )
					{
						var aSplit = sAttr.split('.');
						nNewNode.id = aSplit[0].substr(1, aSplit[0].length-1);
						nNewNode.className = aSplit[1];
					}
					else if ( sAttr.charAt(0) == "#" )
					{
						nNewNode.id = sAttr.substr(1, sAttr.length-1);
					}
					else
					{
						nNewNode.className = sAttr;
					}
	
					i += j; /* Move along the position array */
				}
	
				insert.append( nNewNode );
				insert = $(nNewNode);
			}
			else if ( cOption == '>' )
			{
				/* End container div */
				insert = insert.parent();
			}
			// @todo Move options into their own plugins?
			else if ( cOption == 'l' && features.bPaginate && features.bLengthChange )
			{
				/* Length */
				featureNode = _fnFeatureHtmlLength( oSettings );
			}
			else if ( cOption == 'f' && features.bFilter )
			{
				/* Filter */
				featureNode = _fnFeatureHtmlFilter( oSettings );
			}
			else if ( cOption == 'r' && features.bProcessing )
			{
				/* pRocessing */
				featureNode = _fnFeatureHtmlProcessing( oSettings );
			}
			else if ( cOption == 't' )
			{
				/* Table */
				featureNode = _fnFeatureHtmlTable( oSettings );
			}
			else if ( cOption ==  'i' && features.bInfo )
			{
				/* Info */
				featureNode = _fnFeatureHtmlInfo( oSettings );
			}
			else if ( cOption == 'p' && features.bPaginate )
			{
				/* Pagination */
				featureNode = _fnFeatureHtmlPaginate( oSettings );
			}
			else if ( DataTable.ext.feature.length !== 0 )
			{
				/* Plug-in features */
				var aoFeatures = DataTable.ext.feature;
				for ( var k=0, kLen=aoFeatures.length ; k<kLen ; k++ )
				{
					if ( cOption == aoFeatures[k].cFeature )
					{
						featureNode = aoFeatures[k].fnInit( oSettings );
						break;
					}
				}
			}
	
			/* Add to the 2D features array */
			if ( featureNode )
			{
				var aanFeatures = oSettings.aanFeatures;
	
				if ( ! aanFeatures[cOption] )
				{
					aanFeatures[cOption] = [];
				}
	
				aanFeatures[cOption].push( featureNode );
				insert.append( featureNode );
			}
		}
	
		/* Built our DOM structure - replace the holding div with what we want */
		holding.replaceWith( insert );
		oSettings.nHolding = null;
	}
	
	
	/**
	 * Use the DOM source to create up an array of header cells. The idea here is to
	 * create a layout grid (array) of rows x columns, which contains a reference
	 * to the cell that that point in the grid (regardless of col/rowspan), such that
	 * any column / row could be removed and the new grid constructed
	 *  @param array {object} aLayout Array to store the calculated layout in
	 *  @param {node} nThead The header/footer element for the table
	 *  @memberof DataTable#oApi
	 */
	function _fnDetectHeader ( aLayout, nThead )
	{
		var nTrs = $(nThead).children('tr');
		var nTr, nCell;
		var i, k, l, iLen, jLen, iColShifted, iColumn, iColspan, iRowspan;
		var bUnique;
		var fnShiftCol = function ( a, i, j ) {
			var k = a[i];
	                while ( k[j] ) {
				j++;
			}
			return j;
		};
	
		aLayout.splice( 0, aLayout.length );
	
		/* We know how many rows there are in the layout - so prep it */
		for ( i=0, iLen=nTrs.length ; i<iLen ; i++ )
		{
			aLayout.push( [] );
		}
	
		/* Calculate a layout array */
		for ( i=0, iLen=nTrs.length ; i<iLen ; i++ )
		{
			nTr = nTrs[i];
			iColumn = 0;
	
			/* For every cell in the row... */
			nCell = nTr.firstChild;
			while ( nCell ) {
				if ( nCell.nodeName.toUpperCase() == "TD" ||
				     nCell.nodeName.toUpperCase() == "TH" )
				{
					/* Get the col and rowspan attributes from the DOM and sanitise them */
					iColspan = nCell.getAttribute('colspan') * 1;
					iRowspan = nCell.getAttribute('rowspan') * 1;
					iColspan = (!iColspan || iColspan===0 || iColspan===1) ? 1 : iColspan;
					iRowspan = (!iRowspan || iRowspan===0 || iRowspan===1) ? 1 : iRowspan;
	
					/* There might be colspan cells already in this row, so shift our target
					 * accordingly
					 */
					iColShifted = fnShiftCol( aLayout, i, iColumn );
	
					/* Cache calculation for unique columns */
					bUnique = iColspan === 1 ? true : false;
	
					/* If there is col / rowspan, copy the information into the layout grid */
					for ( l=0 ; l<iColspan ; l++ )
					{
						for ( k=0 ; k<iRowspan ; k++ )
						{
							aLayout[i+k][iColShifted+l] = {
								"cell": nCell,
								"unique": bUnique
							};
							aLayout[i+k].nTr = nTr;
						}
					}
				}
				nCell = nCell.nextSibling;
			}
		}
	}
	
	
	/**
	 * Get an array of unique th elements, one for each column
	 *  @param {object} oSettings dataTables settings object
	 *  @param {node} nHeader automatically detect the layout from this node - optional
	 *  @param {array} aLayout thead/tfoot layout from _fnDetectHeader - optional
	 *  @returns array {node} aReturn list of unique th's
	 *  @memberof DataTable#oApi
	 */
	function _fnGetUniqueThs ( oSettings, nHeader, aLayout )
	{
		var aReturn = [];
		if ( !aLayout )
		{
			aLayout = oSettings.aoHeader;
			if ( nHeader )
			{
				aLayout = [];
				_fnDetectHeader( aLayout, nHeader );
			}
		}
	
		for ( var i=0, iLen=aLayout.length ; i<iLen ; i++ )
		{
			for ( var j=0, jLen=aLayout[i].length ; j<jLen ; j++ )
			{
				if ( aLayout[i][j].unique &&
					 (!aReturn[j] || !oSettings.bSortCellsTop) )
				{
					aReturn[j] = aLayout[i][j].cell;
				}
			}
		}
	
		return aReturn;
	}
	
	/**
	 * Set the start position for draw
	 *  @param {object} oSettings dataTables settings object
	 */
	function _fnStart( oSettings )
	{
		var bServerSide = _fnDataSource( oSettings ) == 'ssp';
		var iInitDisplayStart = oSettings.iInitDisplayStart;
	
		// Check and see if we have an initial draw position from state saving
		if ( iInitDisplayStart !== undefined && iInitDisplayStart !== -1 )
		{
			oSettings._iDisplayStart = bServerSide ?
				iInitDisplayStart :
				iInitDisplayStart >= oSettings.fnRecordsDisplay() ?
					0 :
					iInitDisplayStart;
	
			oSettings.iInitDisplayStart = -1;
		}
	}
	
	/**
	 * Create an Ajax call based on the table's settings, taking into account that
	 * parameters can have multiple forms, and backwards compatibility.
	 *
	 * @param {object} oSettings dataTables settings object
	 * @param {array} data Data to send to the server, required by
	 *     DataTables - may be augmented by developer callbacks
	 * @param {function} fn Callback function to run when data is obtained
	 */
	function _fnBuildAjax( oSettings, data, fn )
	{
		// Compatibility with 1.9-, allow fnServerData and event to manipulate
		_fnCallbackFire( oSettings, 'aoServerParams', 'serverParams', [data] );
	
		// Convert to object based for 1.10+ if using the old array scheme which can
		// come from server-side processing or serverParams
		if ( data && Array.isArray(data) ) {
			var tmp = {};
			var rbracket = /(.*?)\[\]$/;
	
			$.each( data, function (key, val) {
				var match = val.name.match(rbracket);
	
				if ( match ) {
					// Support for arrays
					var name = match[0];
	
					if ( ! tmp[ name ] ) {
						tmp[ name ] = [];
					}
					tmp[ name ].push( val.value );
				}
				else {
					tmp[val.name] = val.value;
				}
			} );
			data = tmp;
		}
	
		var ajaxData;
		var ajax = oSettings.ajax;
		var instance = oSettings.oInstance;
		var callback = function ( json ) {
			var status = oSettings.jqXHR
				? oSettings.jqXHR.status
				: null;
	
			if ( json === null || (typeof status === 'number' && status == 204 ) ) {
				json = {};
				_fnAjaxDataSrc( oSettings, json, [] );
			}
	
			var error = json.error || json.sError;
			if ( error ) {
				_fnLog( oSettings, 0, error );
			}
	
			oSettings.json = json;
	
			_fnCallbackFire( oSettings, null, 'xhr', [oSettings, json, oSettings.jqXHR] );
			fn( json );
		};
	
		if ( $.isPlainObject( ajax ) && ajax.data )
		{
			ajaxData = ajax.data;
	
			var newData = typeof ajaxData === 'function' ?
				ajaxData( data, oSettings ) :  // fn can manipulate data or return
				ajaxData;                      // an object object or array to merge
	
			// If the function returned something, use that alone
			data = typeof ajaxData === 'function' && newData ?
				newData :
				$.extend( true, data, newData );
	
			// Remove the data property as we've resolved it already and don't want
			// jQuery to do it again (it is restored at the end of the function)
			delete ajax.data;
		}
	
		var baseAjax = {
			"data": data,
			"success": callback,
			"dataType": "json",
			"cache": false,
			"type": oSettings.sServerMethod,
			"error": function (xhr, error, thrown) {
				var ret = _fnCallbackFire( oSettings, null, 'xhr', [oSettings, null, oSettings.jqXHR] );
	
				if ( $.inArray( true, ret ) === -1 ) {
					if ( error == "parsererror" ) {
						_fnLog( oSettings, 0, 'Invalid JSON response', 1 );
					}
					else if ( xhr.readyState === 4 ) {
						_fnLog( oSettings, 0, 'Ajax error', 7 );
					}
				}
	
				_fnProcessingDisplay( oSettings, false );
			}
		};
	
		// Store the data submitted for the API
		oSettings.oAjaxData = data;
	
		// Allow plug-ins and external processes to modify the data
		_fnCallbackFire( oSettings, null, 'preXhr', [oSettings, data] );
	
		if ( oSettings.fnServerData )
		{
			// DataTables 1.9- compatibility
			oSettings.fnServerData.call( instance,
				oSettings.sAjaxSource,
				$.map( data, function (val, key) { // Need to convert back to 1.9 trad format
					return { name: key, value: val };
				} ),
				callback,
				oSettings
			);
		}
		else if ( oSettings.sAjaxSource || typeof ajax === 'string' )
		{
			// DataTables 1.9- compatibility
			oSettings.jqXHR = $.ajax( $.extend( baseAjax, {
				url: ajax || oSettings.sAjaxSource
			} ) );
		}
		else if ( typeof ajax === 'function' )
		{
			// Is a function - let the caller define what needs to be done
			oSettings.jqXHR = ajax.call( instance, data, callback, oSettings );
		}
		else
		{
			// Object to extend the base settings
			oSettings.jqXHR = $.ajax( $.extend( baseAjax, ajax ) );
	
			// Restore for next time around
			ajax.data = ajaxData;
		}
	}
	
	
	/**
	 * Update the table using an Ajax call
	 *  @param {object} settings dataTables settings object
	 *  @returns {boolean} Block the table drawing or not
	 *  @memberof DataTable#oApi
	 */
	function _fnAjaxUpdate( settings )
	{
		settings.iDraw++;
		_fnProcessingDisplay( settings, true );
	
		_fnBuildAjax(
			settings,
			_fnAjaxParameters( settings ),
			function(json) {
				_fnAjaxUpdateDraw( settings, json );
			}
		);
	}
	
	
	/**
	 * Build up the parameters in an object needed for a server-side processing
	 * request. Note that this is basically done twice, is different ways - a modern
	 * method which is used by default in DataTables 1.10 which uses objects and
	 * arrays, or the 1.9- method with is name / value pairs. 1.9 method is used if
	 * the sAjaxSource option is used in the initialisation, or the legacyAjax
	 * option is set.
	 *  @param {object} oSettings dataTables settings object
	 *  @returns {bool} block the table drawing or not
	 *  @memberof DataTable#oApi
	 */
	function _fnAjaxParameters( settings )
	{
		var
			columns = settings.aoColumns,
			columnCount = columns.length,
			features = settings.oFeatures,
			preSearch = settings.oPreviousSearch,
			preColSearch = settings.aoPreSearchCols,
			i, data = [], dataProp, column, columnSearch,
			sort = _fnSortFlatten( settings ),
			displayStart = settings._iDisplayStart,
			displayLength = features.bPaginate !== false ?
				settings._iDisplayLength :
				-1;
	
		var param = function ( name, value ) {
			data.push( { 'name': name, 'value': value } );
		};
	
		// DataTables 1.9- compatible method
		param( 'sEcho',          settings.iDraw );
		param( 'iColumns',       columnCount );
		param( 'sColumns',       _pluck( columns, 'sName' ).join(',') );
		param( 'iDisplayStart',  displayStart );
		param( 'iDisplayLength', displayLength );
	
		// DataTables 1.10+ method
		var d = {
			draw:    settings.iDraw,
			columns: [],
			order:   [],
			start:   displayStart,
			length:  displayLength,
			search:  {
				value: preSearch.sSearch,
				regex: preSearch.bRegex
			}
		};
	
		for ( i=0 ; i<columnCount ; i++ ) {
			column = columns[i];
			columnSearch = preColSearch[i];
			dataProp = typeof column.mData=="function" ? 'function' : column.mData ;
	
			d.columns.push( {
				data:       dataProp,
				name:       column.sName,
				searchable: column.bSearchable,
				orderable:  column.bSortable,
				search:     {
					value: columnSearch.sSearch,
					regex: columnSearch.bRegex
				}
			} );
	
			param( "mDataProp_"+i, dataProp );
	
			if ( features.bFilter ) {
				param( 'sSearch_'+i,     columnSearch.sSearch );
				param( 'bRegex_'+i,      columnSearch.bRegex );
				param( 'bSearchable_'+i, column.bSearchable );
			}
	
			if ( features.bSort ) {
				param( 'bSortable_'+i, column.bSortable );
			}
		}
	
		if ( features.bFilter ) {
			param( 'sSearch', preSearch.sSearch );
			param( 'bRegex', preSearch.bRegex );
		}
	
		if ( features.bSort ) {
			$.each( sort, function ( i, val ) {
				d.order.push( { column: val.col, dir: val.dir } );
	
				param( 'iSortCol_'+i, val.col );
				param( 'sSortDir_'+i, val.dir );
			} );
	
			param( 'iSortingCols', sort.length );
		}
	
		// If the legacy.ajax parameter is null, then we automatically decide which
		// form to use, based on sAjaxSource
		var legacy = DataTable.ext.legacy.ajax;
		if ( legacy === null ) {
			return settings.sAjaxSource ? data : d;
		}
	
		// Otherwise, if legacy has been specified then we use that to decide on the
		// form
		return legacy ? data : d;
	}
	
	
	/**
	 * Data the data from the server (nuking the old) and redraw the table
	 *  @param {object} oSettings dataTables settings object
	 *  @param {object} json json data return from the server.
	 *  @param {string} json.sEcho Tracking flag for DataTables to match requests
	 *  @param {int} json.iTotalRecords Number of records in the data set, not accounting for filtering
	 *  @param {int} json.iTotalDisplayRecords Number of records in the data set, accounting for filtering
	 *  @param {array} json.aaData The data to display on this page
	 *  @param {string} [json.sColumns] Column ordering (sName, comma separated)
	 *  @memberof DataTable#oApi
	 */
	function _fnAjaxUpdateDraw ( settings, json )
	{
		// v1.10 uses camelCase variables, while 1.9 uses Hungarian notation.
		// Support both
		var compat = function ( old, modern ) {
			return json[old] !== undefined ? json[old] : json[modern];
		};
	
		var data = _fnAjaxDataSrc( settings, json );
		var draw            = compat( 'sEcho',                'draw' );
		var recordsTotal    = compat( 'iTotalRecords',        'recordsTotal' );
		var recordsFiltered = compat( 'iTotalDisplayRecords', 'recordsFiltered' );
	
		if ( draw !== undefined ) {
			// Protect against out of sequence returns
			if ( draw*1 < settings.iDraw ) {
				return;
			}
			settings.iDraw = draw * 1;
		}
	
		// No data in returned object, so rather than an array, we show an empty table
		if ( ! data ) {
			data = [];
		}
	
		_fnClearTable( settings );
		settings._iRecordsTotal   = parseInt(recordsTotal, 10);
		settings._iRecordsDisplay = parseInt(recordsFiltered, 10);
	
		for ( var i=0, ien=data.length ; i<ien ; i++ ) {
			_fnAddData( settings, data[i] );
		}
		settings.aiDisplay = settings.aiDisplayMaster.slice();
	
		_fnDraw( settings, true );
	
		if ( ! settings._bInitComplete ) {
			_fnInitComplete( settings, json );
		}
	
		_fnProcessingDisplay( settings, false );
	}
	
	
	/**
	 * Get the data from the JSON data source to use for drawing a table. Using
	 * `_fnGetObjectDataFn` allows the data to be sourced from a property of the
	 * source object, or from a processing function.
	 *  @param {object} oSettings dataTables settings object
	 *  @param  {object} json Data source object / array from the server
	 *  @return {array} Array of data to use
	 */
	 function _fnAjaxDataSrc ( oSettings, json, write )
	 {
		var dataSrc = $.isPlainObject( oSettings.ajax ) && oSettings.ajax.dataSrc !== undefined ?
			oSettings.ajax.dataSrc :
			oSettings.sAjaxDataProp; // Compatibility with 1.9-.
	
		if ( ! write ) {
			if ( dataSrc === 'data' ) {
				// If the default, then we still want to support the old style, and safely ignore
				// it if possible
				return json.aaData || json[dataSrc];
			}
	
			return dataSrc !== "" ?
				_fnGetObjectDataFn( dataSrc )( json ) :
				json;
		}
	
		// set
		_fnSetObjectDataFn( dataSrc )( json, write );
	}
	
	/**
	 * Generate the node required for filtering text
	 *  @returns {node} Filter control element
	 *  @param {object} oSettings dataTables settings object
	 *  @memberof DataTable#oApi
	 */
	function _fnFeatureHtmlFilter ( settings )
	{
		var classes = settings.oClasses;
		var tableId = settings.sTableId;
		var language = settings.oLanguage;
		var previousSearch = settings.oPreviousSearch;
		var features = settings.aanFeatures;
		var input = '<input type="search" class="'+classes.sFilterInput+'"/>';
	
		var str = language.sSearch;
		str = str.match(/_INPUT_/) ?
			str.replace('_INPUT_', input) :
			str+input;
	
		var filter = $('<div/>', {
				'id': ! features.f ? tableId+'_filter' : null,
				'class': classes.sFilter
			} )
			.append( $('<label/>' ).append( str ) );
	
		var searchFn = function(event) {
			/* Update all other filter input elements for the new display */
			var n = features.f;
			var val = !this.value ? "" : this.value; // mental IE8 fix :-(
			if(previousSearch.return && event.key !== "Enter") {
				return;
			}
			/* Now do the filter */
			if ( val != previousSearch.sSearch ) {
				_fnFilterComplete( settings, {
					"sSearch": val,
					"bRegex": previousSearch.bRegex,
					"bSmart": previousSearch.bSmart ,
					"bCaseInsensitive": previousSearch.bCaseInsensitive,
					"return": previousSearch.return
				} );
	
				// Need to redraw, without resorting
				settings._iDisplayStart = 0;
				_fnDraw( settings );
			}
		};
	
		var searchDelay = settings.searchDelay !== null ?
			settings.searchDelay :
			_fnDataSource( settings ) === 'ssp' ?
				400 :
				0;
	
		var jqFilter = $('input', filter)
			.val( previousSearch.sSearch )
			.attr( 'placeholder', language.sSearchPlaceholder )
			.on(
				'keyup.DT search.DT input.DT paste.DT cut.DT',
				searchDelay ?
					_fnThrottle( searchFn, searchDelay ) :
					searchFn
			)
			.on( 'mouseup', function(e) {
				// Edge fix! Edge 17 does not trigger anything other than mouse events when clicking
				// on the clear icon (Edge bug 17584515). This is safe in other browsers as `searchFn`
				// checks the value to see if it has changed. In other browsers it won't have.
				setTimeout( function () {
					searchFn.call(jqFilter[0], e);
				}, 10);
			} )
			.on( 'keypress.DT', function(e) {
				/* Prevent form submission */
				if ( e.keyCode == 13 ) {
					return false;
				}
			} )
			.attr('aria-controls', tableId);
	
		// Update the input elements whenever the table is filtered
		$(settings.nTable).on( 'search.dt.DT', function ( ev, s ) {
			if ( settings === s ) {
				// IE9 throws an 'unknown error' if document.activeElement is used
				// inside an iframe or frame...
				try {
					if ( jqFilter[0] !== document.activeElement ) {
						jqFilter.val( previousSearch.sSearch );
					}
				}
				catch ( e ) {}
			}
		} );
	
		return filter[0];
	}
	
	
	/**
	 * Filter the table using both the global filter and column based filtering
	 *  @param {object} oSettings dataTables settings object
	 *  @param {object} oSearch search information
	 *  @param {int} [iForce] force a research of the master array (1) or not (undefined or 0)
	 *  @memberof DataTable#oApi
	 */
	function _fnFilterComplete ( oSettings, oInput, iForce )
	{
		var oPrevSearch = oSettings.oPreviousSearch;
		var aoPrevSearch = oSettings.aoPreSearchCols;
		var fnSaveFilter = function ( oFilter ) {
			/* Save the filtering values */
			oPrevSearch.sSearch = oFilter.sSearch;
			oPrevSearch.bRegex = oFilter.bRegex;
			oPrevSearch.bSmart = oFilter.bSmart;
			oPrevSearch.bCaseInsensitive = oFilter.bCaseInsensitive;
			oPrevSearch.return = oFilter.return;
		};
		var fnRegex = function ( o ) {
			// Backwards compatibility with the bEscapeRegex option
			return o.bEscapeRegex !== undefined ? !o.bEscapeRegex : o.bRegex;
		};
	
		// Resolve any column types that are unknown due to addition or invalidation
		// @todo As per sort - can this be moved into an event handler?
		_fnColumnTypes( oSettings );
	
		/* In server-side processing all filtering is done by the server, so no point hanging around here */
		if ( _fnDataSource( oSettings ) != 'ssp' )
		{
			/* Global filter */
			_fnFilter( oSettings, oInput.sSearch, iForce, fnRegex(oInput), oInput.bSmart, oInput.bCaseInsensitive, oInput.return );
			fnSaveFilter( oInput );
	
			/* Now do the individual column filter */
			for ( var i=0 ; i<aoPrevSearch.length ; i++ )
			{
				_fnFilterColumn( oSettings, aoPrevSearch[i].sSearch, i, fnRegex(aoPrevSearch[i]),
					aoPrevSearch[i].bSmart, aoPrevSearch[i].bCaseInsensitive );
			}
	
			/* Custom filtering */
			_fnFilterCustom( oSettings );
		}
		else
		{
			fnSaveFilter( oInput );
		}
	
		/* Tell the draw function we have been filtering */
		oSettings.bFiltered = true;
		_fnCallbackFire( oSettings, null, 'search', [oSettings] );
	}
	
	
	/**
	 * Apply custom filtering functions
	 *  @param {object} oSettings dataTables settings object
	 *  @memberof DataTable#oApi
	 */
	function _fnFilterCustom( settings )
	{
		var filters = DataTable.ext.search;
		var displayRows = settings.aiDisplay;
		var row, rowIdx;
	
		for ( var i=0, ien=filters.length ; i<ien ; i++ ) {
			var rows = [];
	
			// Loop over each row and see if it should be included
			for ( var j=0, jen=displayRows.length ; j<jen ; j++ ) {
				rowIdx = displayRows[ j ];
				row = settings.aoData[ rowIdx ];
	
				if ( filters[i]( settings, row._aFilterData, rowIdx, row._aData, j ) ) {
					rows.push( rowIdx );
				}
			}
	
			// So the array reference doesn't break set the results into the
			// existing array
			displayRows.length = 0;
			$.merge( displayRows, rows );
		}
	}
	
	
	/**
	 * Filter the table on a per-column basis
	 *  @param {object} oSettings dataTables settings object
	 *  @param {string} sInput string to filter on
	 *  @param {int} iColumn column to filter
	 *  @param {bool} bRegex treat search string as a regular expression or not
	 *  @param {bool} bSmart use smart filtering or not
	 *  @param {bool} bCaseInsensitive Do case insensitive matching or not
	 *  @memberof DataTable#oApi
	 */
	function _fnFilterColumn ( settings, searchStr, colIdx, regex, smart, caseInsensitive )
	{
		if ( searchStr === '' ) {
			return;
		}
	
		var data;
		var out = [];
		var display = settings.aiDisplay;
		var rpSearch = _fnFilterCreateSearch( searchStr, regex, smart, caseInsensitive );
	
		for ( var i=0 ; i<display.length ; i++ ) {
			data = settings.aoData[ display[i] ]._aFilterData[ colIdx ];
	
			if ( rpSearch.test( data ) ) {
				out.push( display[i] );
			}
		}
	
		settings.aiDisplay = out;
	}
	
	
	/**
	 * Filter the data table based on user input and draw the table
	 *  @param {object} settings dataTables settings object
	 *  @param {string} input string to filter on
	 *  @param {int} force optional - force a research of the master array (1) or not (undefined or 0)
	 *  @param {bool} regex treat as a regular expression or not
	 *  @param {bool} smart perform smart filtering or not
	 *  @param {bool} caseInsensitive Do case insensitive matching or not
	 *  @memberof DataTable#oApi
	 */
	function _fnFilter( settings, input, force, regex, smart, caseInsensitive )
	{
		var rpSearch = _fnFilterCreateSearch( input, regex, smart, caseInsensitive );
		var prevSearch = settings.oPreviousSearch.sSearch;
		var displayMaster = settings.aiDisplayMaster;
		var display, invalidated, i;
		var filtered = [];
	
		// Need to take account of custom filtering functions - always filter
		if ( DataTable.ext.search.length !== 0 ) {
			force = true;
		}
	
		// Check if any of the rows were invalidated
		invalidated = _fnFilterData( settings );
	
		// If the input is blank - we just want the full data set
		if ( input.length <= 0 ) {
			settings.aiDisplay = displayMaster.slice();
		}
		else {
			// New search - start from the master array
			if ( invalidated ||
				 force ||
				 regex ||
				 prevSearch.length > input.length ||
				 input.indexOf(prevSearch) !== 0 ||
				 settings.bSorted // On resort, the display master needs to be
				                  // re-filtered since indexes will have changed
			) {
				settings.aiDisplay = displayMaster.slice();
			}
	
			// Search the display array
			display = settings.aiDisplay;
	
			for ( i=0 ; i<display.length ; i++ ) {
				if ( rpSearch.test( settings.aoData[ display[i] ]._sFilterRow ) ) {
					filtered.push( display[i] );
				}
			}
	
			settings.aiDisplay = filtered;
		}
	}
	
	
	/**
	 * Build a regular expression object suitable for searching a table
	 *  @param {string} sSearch string to search for
	 *  @param {bool} bRegex treat as a regular expression or not
	 *  @param {bool} bSmart perform smart filtering or not
	 *  @param {bool} bCaseInsensitive Do case insensitive matching or not
	 *  @returns {RegExp} constructed object
	 *  @memberof DataTable#oApi
	 */
	function _fnFilterCreateSearch( search, regex, smart, caseInsensitive )
	{
		search = regex ?
			search :
			_fnEscapeRegex( search );
		
		if ( smart ) {
			/* For smart filtering we want to allow the search to work regardless of
			 * word order. We also want double quoted text to be preserved, so word
			 * order is important - a la google. So this is what we want to
			 * generate:
			 * 
			 * ^(?=.*?\bone\b)(?=.*?\btwo three\b)(?=.*?\bfour\b).*$
			 */
			var a = $.map( search.match( /"[^"]+"|[^ ]+/g ) || [''], function ( word ) {
				if ( word.charAt(0) === '"' ) {
					var m = word.match( /^"(.*)"$/ );
					word = m ? m[1] : word;
				}
	
				return word.replace('"', '');
			} );
	
			search = '^(?=.*?'+a.join( ')(?=.*?' )+').*$';
		}
	
		return new RegExp( search, caseInsensitive ? 'i' : '' );
	}
	
	
	/**
	 * Escape a string such that it can be used in a regular expression
	 *  @param {string} sVal string to escape
	 *  @returns {string} escaped string
	 *  @memberof DataTable#oApi
	 */
	var _fnEscapeRegex = DataTable.util.escapeRegex;
	
	var __filter_div = $('<div>')[0];
	var __filter_div_textContent = __filter_div.textContent !== undefined;
	
	// Update the filtering data for each row if needed (by invalidation or first run)
	function _fnFilterData ( settings )
	{
		var columns = settings.aoColumns;
		var column;
		var i, j, ien, jen, filterData, cellData, row;
		var wasInvalidated = false;
	
		for ( i=0, ien=settings.aoData.length ; i<ien ; i++ ) {
			row = settings.aoData[i];
	
			if ( ! row._aFilterData ) {
				filterData = [];
	
				for ( j=0, jen=columns.length ; j<jen ; j++ ) {
					column = columns[j];
	
					if ( column.bSearchable ) {
						cellData = _fnGetCellData( settings, i, j, 'filter' );
	
						// Search in DataTables 1.10 is string based. In 1.11 this
						// should be altered to also allow strict type checking.
						if ( cellData === null ) {
							cellData = '';
						}
	
						if ( typeof cellData !== 'string' && cellData.toString ) {
							cellData = cellData.toString();
						}
					}
					else {
						cellData = '';
					}
	
					// If it looks like there is an HTML entity in the string,
					// attempt to decode it so sorting works as expected. Note that
					// we could use a single line of jQuery to do this, but the DOM
					// method used here is much faster http://jsperf.com/html-decode
					if ( cellData.indexOf && cellData.indexOf('&') !== -1 ) {
						__filter_div.innerHTML = cellData;
						cellData = __filter_div_textContent ?
							__filter_div.textContent :
							__filter_div.innerText;
					}
	
					if ( cellData.replace ) {
						cellData = cellData.replace(/[\r\n\u2028]/g, '');
					}
	
					filterData.push( cellData );
				}
	
				row._aFilterData = filterData;
				row._sFilterRow = filterData.join('  ');
				wasInvalidated = true;
			}
		}
	
		return wasInvalidated;
	}
	
	
	/**
	 * Convert from the internal Hungarian notation to camelCase for external
	 * interaction
	 *  @param {object} obj Object to convert
	 *  @returns {object} Inverted object
	 *  @memberof DataTable#oApi
	 */
	function _fnSearchToCamel ( obj )
	{
		return {
			search:          obj.sSearch,
			smart:           obj.bSmart,
			regex:           obj.bRegex,
			caseInsensitive: obj.bCaseInsensitive
		};
	}
	
	
	
	/**
	 * Convert from camelCase notation to the internal Hungarian. We could use the
	 * Hungarian convert function here, but this is cleaner
	 *  @param {object} obj Object to convert
	 *  @returns {object} Inverted object
	 *  @memberof DataTable#oApi
	 */
	function _fnSearchToHung ( obj )
	{
		return {
			sSearch:          obj.search,
			bSmart:           obj.smart,
			bRegex:           obj.regex,
			bCaseInsensitive: obj.caseInsensitive
		};
	}
	
	/**
	 * Generate the node required for the info display
	 *  @param {object} oSettings dataTables settings object
	 *  @returns {node} Information element
	 *  @memberof DataTable#oApi
	 */
	function _fnFeatureHtmlInfo ( settings )
	{
		var
			tid = settings.sTableId,
			nodes = settings.aanFeatures.i,
			n = $('<div/>', {
				'class': settings.oClasses.sInfo,
				'id': ! nodes ? tid+'_info' : null
			} );
	
		if ( ! nodes ) {
			// Update display on each draw
			settings.aoDrawCallback.push( {
				"fn": _fnUpdateInfo,
				"sName": "information"
			} );
	
			n
				.attr( 'role', 'status' )
				.attr( 'aria-live', 'polite' );
	
			// Table is described by our info div
			$(settings.nTable).attr( 'aria-describedby', tid+'_info' );
		}
	
		return n[0];
	}
	
	
	/**
	 * Update the information elements in the display
	 *  @param {object} settings dataTables settings object
	 *  @memberof DataTable#oApi
	 */
	function _fnUpdateInfo ( settings )
	{
		/* Show information about the table */
		var nodes = settings.aanFeatures.i;
		if ( nodes.length === 0 ) {
			return;
		}
	
		var
			lang  = settings.oLanguage,
			start = settings._iDisplayStart+1,
			end   = settings.fnDisplayEnd(),
			max   = settings.fnRecordsTotal(),
			total = settings.fnRecordsDisplay(),
			out   = total ?
				lang.sInfo :
				lang.sInfoEmpty;
	
		if ( total !== max ) {
			/* Record set after filtering */
			out += ' ' + lang.sInfoFiltered;
		}
	
		// Convert the macros
		out += lang.sInfoPostFix;
		out = _fnInfoMacros( settings, out );
	
		var callback = lang.fnInfoCallback;
		if ( callback !== null ) {
			out = callback.call( settings.oInstance,
				settings, start, end, max, total, out
			);
		}
	
		$(nodes).html( out );
	}
	
	
	function _fnInfoMacros ( settings, str )
	{
		// When infinite scrolling, we are always starting at 1. _iDisplayStart is used only
		// internally
		var
			formatter  = settings.fnFormatNumber,
			start      = settings._iDisplayStart+1,
			len        = settings._iDisplayLength,
			vis        = settings.fnRecordsDisplay(),
			all        = len === -1;
	
		return str.
			replace(/_START_/g, formatter.call( settings, start ) ).
			replace(/_END_/g,   formatter.call( settings, settings.fnDisplayEnd() ) ).
			replace(/_MAX_/g,   formatter.call( settings, settings.fnRecordsTotal() ) ).
			replace(/_TOTAL_/g, formatter.call( settings, vis ) ).
			replace(/_PAGE_/g,  formatter.call( settings, all ? 1 : Math.ceil( start / len ) ) ).
			replace(/_PAGES_/g, formatter.call( settings, all ? 1 : Math.ceil( vis / len ) ) );
	}
	
	
	
	/**
	 * Draw the table for the first time, adding all required features
	 *  @param {object} settings dataTables settings object
	 *  @memberof DataTable#oApi
	 */
	function _fnInitialise ( settings )
	{
		var i, iLen, iAjaxStart=settings.iInitDisplayStart;
		var columns = settings.aoColumns, column;
		var features = settings.oFeatures;
		var deferLoading = settings.bDeferLoading; // value modified by the draw
	
		/* Ensure that the table data is fully initialised */
		if ( ! settings.bInitialised ) {
			setTimeout( function(){ _fnInitialise( settings ); }, 200 );
			return;
		}
	
		/* Show the display HTML options */
		_fnAddOptionsHtml( settings );
	
		/* Build and draw the header / footer for the table */
		_fnBuildHead( settings );
		_fnDrawHead( settings, settings.aoHeader );
		_fnDrawHead( settings, settings.aoFooter );
	
		/* Okay to show that something is going on now */
		_fnProcessingDisplay( settings, true );
	
		/* Calculate sizes for columns */
		if ( features.bAutoWidth ) {
			_fnCalculateColumnWidths( settings );
		}
	
		for ( i=0, iLen=columns.length ; i<iLen ; i++ ) {
			column = columns[i];
	
			if ( column.sWidth ) {
				column.nTh.style.width = _fnStringToCss( column.sWidth );
			}
		}
	
		_fnCallbackFire( settings, null, 'preInit', [settings] );
	
		// If there is default sorting required - let's do it. The sort function
		// will do the drawing for us. Otherwise we draw the table regardless of the
		// Ajax source - this allows the table to look initialised for Ajax sourcing
		// data (show 'loading' message possibly)
		_fnReDraw( settings );
	
		// Server-side processing init complete is done by _fnAjaxUpdateDraw
		var dataSrc = _fnDataSource( settings );
		if ( dataSrc != 'ssp' || deferLoading ) {
			// if there is an ajax source load the data
			if ( dataSrc == 'ajax' ) {
				_fnBuildAjax( settings, [], function(json) {
					var aData = _fnAjaxDataSrc( settings, json );
	
					// Got the data - add it to the table
					for ( i=0 ; i<aData.length ; i++ ) {
						_fnAddData( settings, aData[i] );
					}
	
					// Reset the init display for cookie saving. We've already done
					// a filter, and therefore cleared it before. So we need to make
					// it appear 'fresh'
					settings.iInitDisplayStart = iAjaxStart;
	
					_fnReDraw( settings );
	
					_fnProcessingDisplay( settings, false );
					_fnInitComplete( settings, json );
				}, settings );
			}
			else {
				_fnProcessingDisplay( settings, false );
				_fnInitComplete( settings );
			}
		}
	}
	
	
	/**
	 * Draw the table for the first time, adding all required features
	 *  @param {object} oSettings dataTables settings object
	 *  @param {object} [json] JSON from the server that completed the table, if using Ajax source
	 *    with client-side processing (optional)
	 *  @memberof DataTable#oApi
	 */
	function _fnInitComplete ( settings, json )
	{
		settings._bInitComplete = true;
	
		// When data was added after the initialisation (data or Ajax) we need to
		// calculate the column sizing
		if ( json || settings.oInit.aaData ) {
			_fnAdjustColumnSizing( settings );
		}
	
		_fnCallbackFire( settings, null, 'plugin-init', [settings, json] );
		_fnCallbackFire( settings, 'aoInitComplete', 'init', [settings, json] );
	}
	
	
	function _fnLengthChange ( settings, val )
	{
		var len = parseInt( val, 10 );
		settings._iDisplayLength = len;
	
		_fnLengthOverflow( settings );
	
		// Fire length change event
		_fnCallbackFire( settings, null, 'length', [settings, len] );
	}
	
	
	/**
	 * Generate the node required for user display length changing
	 *  @param {object} settings dataTables settings object
	 *  @returns {node} Display length feature node
	 *  @memberof DataTable#oApi
	 */
	function _fnFeatureHtmlLength ( settings )
	{
		var
			classes  = settings.oClasses,
			tableId  = settings.sTableId,
			menu     = settings.aLengthMenu,
			d2       = Array.isArray( menu[0] ),
			lengths  = d2 ? menu[0] : menu,
			language = d2 ? menu[1] : menu;
	
		var select = $('<select/>', {
			'name':          tableId+'_length',
			'aria-controls': tableId,
			'class':         classes.sLengthSelect
		} );
	
		for ( var i=0, ien=lengths.length ; i<ien ; i++ ) {
			select[0][ i ] = new Option(
				typeof language[i] === 'number' ?
					settings.fnFormatNumber( language[i] ) :
					language[i],
				lengths[i]
			);
		}
	
		var div = $('<div><label/></div>').addClass( classes.sLength );
		if ( ! settings.aanFeatures.l ) {
			div[0].id = tableId+'_length';
		}
	
		div.children().append(
			settings.oLanguage.sLengthMenu.replace( '_MENU_', select[0].outerHTML )
		);
	
		// Can't use `select` variable as user might provide their own and the
		// reference is broken by the use of outerHTML
		$('select', div)
			.val( settings._iDisplayLength )
			.on( 'change.DT', function(e) {
				_fnLengthChange( settings, $(this).val() );
				_fnDraw( settings );
			} );
	
		// Update node value whenever anything changes the table's length
		$(settings.nTable).on( 'length.dt.DT', function (e, s, len) {
			if ( settings === s ) {
				$('select', div).val( len );
			}
		} );
	
		return div[0];
	}
	
	
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Note that most of the paging logic is done in
	 * DataTable.ext.pager
	 */
	
	/**
	 * Generate the node required for default pagination
	 *  @param {object} oSettings dataTables settings object
	 *  @returns {node} Pagination feature node
	 *  @memberof DataTable#oApi
	 */
	function _fnFeatureHtmlPaginate ( settings )
	{
		var
			type   = settings.sPaginationType,
			plugin = DataTable.ext.pager[ type ],
			modern = typeof plugin === 'function',
			redraw = function( settings ) {
				_fnDraw( settings );
			},
			node = $('<div/>').addClass( settings.oClasses.sPaging + type )[0],
			features = settings.aanFeatures;
	
		if ( ! modern ) {
			plugin.fnInit( settings, node, redraw );
		}
	
		/* Add a draw callback for the pagination on first instance, to update the paging display */
		if ( ! features.p )
		{
			node.id = settings.sTableId+'_paginate';
	
			settings.aoDrawCallback.push( {
				"fn": function( settings ) {
					if ( modern ) {
						var
							start      = settings._iDisplayStart,
							len        = settings._iDisplayLength,
							visRecords = settings.fnRecordsDisplay(),
							all        = len === -1,
							page = all ? 0 : Math.ceil( start / len ),
							pages = all ? 1 : Math.ceil( visRecords / len ),
							buttons = plugin(page, pages),
							i, ien;
	
						for ( i=0, ien=features.p.length ; i<ien ; i++ ) {
							_fnRenderer( settings, 'pageButton' )(
								settings, features.p[i], i, buttons, page, pages
							);
						}
					}
					else {
						plugin.fnUpdate( settings, redraw );
					}
				},
				"sName": "pagination"
			} );
		}
	
		return node;
	}
	
	
	/**
	 * Alter the display settings to change the page
	 *  @param {object} settings DataTables settings object
	 *  @param {string|int} action Paging action to take: "first", "previous",
	 *    "next" or "last" or page number to jump to (integer)
	 *  @param [bool] redraw Automatically draw the update or not
	 *  @returns {bool} true page has changed, false - no change
	 *  @memberof DataTable#oApi
	 */
	function _fnPageChange ( settings, action, redraw )
	{
		var
			start     = settings._iDisplayStart,
			len       = settings._iDisplayLength,
			records   = settings.fnRecordsDisplay();
	
		if ( records === 0 || len === -1 )
		{
			start = 0;
		}
		else if ( typeof action === "number" )
		{
			start = action * len;
	
			if ( start > records )
			{
				start = 0;
			}
		}
		else if ( action == "first" )
		{
			start = 0;
		}
		else if ( action == "previous" )
		{
			start = len >= 0 ?
				start - len :
				0;
	
			if ( start < 0 )
			{
			  start = 0;
			}
		}
		else if ( action == "next" )
		{
			if ( start + len < records )
			{
				start += len;
			}
		}
		else if ( action == "last" )
		{
			start = Math.floor( (records-1) / len) * len;
		}
		else
		{
			_fnLog( settings, 0, "Unknown paging action: "+action, 5 );
		}
	
		var changed = settings._iDisplayStart !== start;
		settings._iDisplayStart = start;
	
		if ( changed ) {
			_fnCallbackFire( settings, null, 'page', [settings] );
	
			if ( redraw ) {
				_fnDraw( settings );
			}
		}
		else {
			// No change event - paging was called, but no change
			_fnCallbackFire( settings, null, 'page-nc', [settings] );
		}
	
		return changed;
	}
	
	
	
	/**
	 * Generate the node required for the processing node
	 *  @param {object} settings dataTables settings object
	 *  @returns {node} Processing element
	 *  @memberof DataTable#oApi
	 */
	function _fnFeatureHtmlProcessing ( settings )
	{
		return $('<div/>', {
				'id': ! settings.aanFeatures.r ? settings.sTableId+'_processing' : null,
				'class': settings.oClasses.sProcessing,
				'role': 'status'
			} )
			.html( settings.oLanguage.sProcessing )
			.append('<div><div></div><div></div><div></div><div></div></div>')
			.insertBefore( settings.nTable )[0];
	}
	
	
	/**
	 * Display or hide the processing indicator
	 *  @param {object} settings dataTables settings object
	 *  @param {bool} show Show the processing indicator (true) or not (false)
	 *  @memberof DataTable#oApi
	 */
	function _fnProcessingDisplay ( settings, show )
	{
		if ( settings.oFeatures.bProcessing ) {
			$(settings.aanFeatures.r).css( 'display', show ? 'block' : 'none' );
		}
	
		_fnCallbackFire( settings, null, 'processing', [settings, show] );
	}
	
	/**
	 * Add any control elements for the table - specifically scrolling
	 *  @param {object} settings dataTables settings object
	 *  @returns {node} Node to add to the DOM
	 *  @memberof DataTable#oApi
	 */
	function _fnFeatureHtmlTable ( settings )
	{
		var table = $(settings.nTable);
	
		// Scrolling from here on in
		var scroll = settings.oScroll;
	
		if ( scroll.sX === '' && scroll.sY === '' ) {
			return settings.nTable;
		}
	
		var scrollX = scroll.sX;
		var scrollY = scroll.sY;
		var classes = settings.oClasses;
		var caption = table.children('caption');
		var captionSide = caption.length ? caption[0]._captionSide : null;
		var headerClone = $( table[0].cloneNode(false) );
		var footerClone = $( table[0].cloneNode(false) );
		var footer = table.children('tfoot');
		var _div = '<div/>';
		var size = function ( s ) {
			return !s ? null : _fnStringToCss( s );
		};
	
		if ( ! footer.length ) {
			footer = null;
		}
	
		/*
		 * The HTML structure that we want to generate in this function is:
		 *  div - scroller
		 *    div - scroll head
		 *      div - scroll head inner
		 *        table - scroll head table
		 *          thead - thead
		 *    div - scroll body
		 *      table - table (master table)
		 *        thead - thead clone for sizing
		 *        tbody - tbody
		 *    div - scroll foot
		 *      div - scroll foot inner
		 *        table - scroll foot table
		 *          tfoot - tfoot
		 */
		var scroller = $( _div, { 'class': classes.sScrollWrapper } )
			.append(
				$(_div, { 'class': classes.sScrollHead } )
					.css( {
						overflow: 'hidden',
						position: 'relative',
						border: 0,
						width: scrollX ? size(scrollX) : '100%'
					} )
					.append(
						$(_div, { 'class': classes.sScrollHeadInner } )
							.css( {
								'box-sizing': 'content-box',
								width: scroll.sXInner || '100%'
							} )
							.append(
								headerClone
									.removeAttr('id')
									.css( 'margin-left', 0 )
									.append( captionSide === 'top' ? caption : null )
									.append(
										table.children('thead')
									)
							)
					)
			)
			.append(
				$(_div, { 'class': classes.sScrollBody } )
					.css( {
						position: 'relative',
						overflow: 'auto',
						width: size( scrollX )
					} )
					.append( table )
			);
	
		if ( footer ) {
			scroller.append(
				$(_div, { 'class': classes.sScrollFoot } )
					.css( {
						overflow: 'hidden',
						border: 0,
						width: scrollX ? size(scrollX) : '100%'
					} )
					.append(
						$(_div, { 'class': classes.sScrollFootInner } )
							.append(
								footerClone
									.removeAttr('id')
									.css( 'margin-left', 0 )
									.append( captionSide === 'bottom' ? caption : null )
									.append(
										table.children('tfoot')
									)
							)
					)
			);
		}
	
		var children = scroller.children();
		var scrollHead = children[0];
		var scrollBody = children[1];
		var scrollFoot = footer ? children[2] : null;
	
		// When the body is scrolled, then we also want to scroll the headers
		if ( scrollX ) {
			$(scrollBody).on( 'scroll.DT', function (e) {
				var scrollLeft = this.scrollLeft;
	
				scrollHead.scrollLeft = scrollLeft;
	
				if ( footer ) {
					scrollFoot.scrollLeft = scrollLeft;
				}
			} );
		}
	
		$(scrollBody).css('max-height', scrollY);
		if (! scroll.bCollapse) {
			$(scrollBody).css('height', scrollY);
		}
	
		settings.nScrollHead = scrollHead;
		settings.nScrollBody = scrollBody;
		settings.nScrollFoot = scrollFoot;
	
		// On redraw - align columns
		settings.aoDrawCallback.push( {
			"fn": _fnScrollDraw,
			"sName": "scrolling"
		} );
	
		return scroller[0];
	}
	
	
	
	/**
	 * Update the header, footer and body tables for resizing - i.e. column
	 * alignment.
	 *
	 * Welcome to the most horrible function DataTables. The process that this
	 * function follows is basically:
	 *   1. Re-create the table inside the scrolling div
	 *   2. Take live measurements from the DOM
	 *   3. Apply the measurements to align the columns
	 *   4. Clean up
	 *
	 *  @param {object} settings dataTables settings object
	 *  @memberof DataTable#oApi
	 */
	function _fnScrollDraw ( settings )
	{
		// Given that this is such a monster function, a lot of variables are use
		// to try and keep the minimised size as small as possible
		var
			scroll         = settings.oScroll,
			scrollX        = scroll.sX,
			scrollXInner   = scroll.sXInner,
			scrollY        = scroll.sY,
			barWidth       = scroll.iBarWidth,
			divHeader      = $(settings.nScrollHead),
			divHeaderStyle = divHeader[0].style,
			divHeaderInner = divHeader.children('div'),
			divHeaderInnerStyle = divHeaderInner[0].style,
			divHeaderTable = divHeaderInner.children('table'),
			divBodyEl      = settings.nScrollBody,
			divBody        = $(divBodyEl),
			divBodyStyle   = divBodyEl.style,
			divFooter      = $(settings.nScrollFoot),
			divFooterInner = divFooter.children('div'),
			divFooterTable = divFooterInner.children('table'),
			header         = $(settings.nTHead),
			table          = $(settings.nTable),
			tableEl        = table[0],
			tableStyle     = tableEl.style,
			footer         = settings.nTFoot ? $(settings.nTFoot) : null,
			browser        = settings.oBrowser,
			ie67           = browser.bScrollOversize,
			dtHeaderCells  = _pluck( settings.aoColumns, 'nTh' ),
			headerTrgEls, footerTrgEls,
			headerSrcEls, footerSrcEls,
			headerCopy, footerCopy,
			headerWidths=[], footerWidths=[],
			headerContent=[], footerContent=[],
			idx, correction, sanityWidth,
			zeroOut = function(nSizer) {
				var style = nSizer.style;
				style.paddingTop = "0";
				style.paddingBottom = "0";
				style.borderTopWidth = "0";
				style.borderBottomWidth = "0";
				style.height = 0;
			};
	
		// If the scrollbar visibility has changed from the last draw, we need to
		// adjust the column sizes as the table width will have changed to account
		// for the scrollbar
		var scrollBarVis = divBodyEl.scrollHeight > divBodyEl.clientHeight;
		
		if ( settings.scrollBarVis !== scrollBarVis && settings.scrollBarVis !== undefined ) {
			settings.scrollBarVis = scrollBarVis;
			_fnAdjustColumnSizing( settings );
			return; // adjust column sizing will call this function again
		}
		else {
			settings.scrollBarVis = scrollBarVis;
		}
	
		/*
		 * 1. Re-create the table inside the scrolling div
		 */
	
		// Remove the old minimised thead and tfoot elements in the inner table
		table.children('thead, tfoot').remove();
	
		if ( footer ) {
			footerCopy = footer.clone().prependTo( table );
			footerTrgEls = footer.find('tr'); // the original tfoot is in its own table and must be sized
			footerSrcEls = footerCopy.find('tr');
			footerCopy.find('[id]').removeAttr('id');
		}
	
		// Clone the current header and footer elements and then place it into the inner table
		headerCopy = header.clone().prependTo( table );
		headerTrgEls = header.find('tr'); // original header is in its own table
		headerSrcEls = headerCopy.find('tr');
		headerCopy.find('th, td').removeAttr('tabindex');
		headerCopy.find('[id]').removeAttr('id');
	
	
		/*
		 * 2. Take live measurements from the DOM - do not alter the DOM itself!
		 */
	
		// Remove old sizing and apply the calculated column widths
		// Get the unique column headers in the newly created (cloned) header. We want to apply the
		// calculated sizes to this header
		if ( ! scrollX )
		{
			divBodyStyle.width = '100%';
			divHeader[0].style.width = '100%';
		}
	
		$.each( _fnGetUniqueThs( settings, headerCopy ), function ( i, el ) {
			idx = _fnVisibleToColumnIndex( settings, i );
			el.style.width = settings.aoColumns[idx].sWidth;
		} );
	
		if ( footer ) {
			_fnApplyToChildren( function(n) {
				n.style.width = "";
			}, footerSrcEls );
		}
	
		// Size the table as a whole
		sanityWidth = table.outerWidth();
		if ( scrollX === "" ) {
			// No x scrolling
			tableStyle.width = "100%";
	
			// IE7 will make the width of the table when 100% include the scrollbar
			// - which is shouldn't. When there is a scrollbar we need to take this
			// into account.
			if ( ie67 && (table.find('tbody').height() > divBodyEl.offsetHeight ||
				divBody.css('overflow-y') == "scroll")
			) {
				tableStyle.width = _fnStringToCss( table.outerWidth() - barWidth);
			}
	
			// Recalculate the sanity width
			sanityWidth = table.outerWidth();
		}
		else if ( scrollXInner !== "" ) {
			// legacy x scroll inner has been given - use it
			tableStyle.width = _fnStringToCss(scrollXInner);
	
			// Recalculate the sanity width
			sanityWidth = table.outerWidth();
		}
	
		// Hidden header should have zero height, so remove padding and borders. Then
		// set the width based on the real headers
	
		// Apply all styles in one pass
		_fnApplyToChildren( zeroOut, headerSrcEls );
	
		// Read all widths in next pass
		_fnApplyToChildren( function(nSizer) {
			var style = window.getComputedStyle ?
				window.getComputedStyle(nSizer).width :
				_fnStringToCss( $(nSizer).width() );
	
			headerContent.push( nSizer.innerHTML );
			headerWidths.push( style );
		}, headerSrcEls );
	
		// Apply all widths in final pass
		_fnApplyToChildren( function(nToSize, i) {
			nToSize.style.width = headerWidths[i];
		}, headerTrgEls );
	
		$(headerSrcEls).css('height', 0);
	
		/* Same again with the footer if we have one */
		if ( footer )
		{
			_fnApplyToChildren( zeroOut, footerSrcEls );
	
			_fnApplyToChildren( function(nSizer) {
				footerContent.push( nSizer.innerHTML );
				footerWidths.push( _fnStringToCss( $(nSizer).css('width') ) );
			}, footerSrcEls );
	
			_fnApplyToChildren( function(nToSize, i) {
				nToSize.style.width = footerWidths[i];
			}, footerTrgEls );
	
			$(footerSrcEls).height(0);
		}
	
	
		/*
		 * 3. Apply the measurements
		 */
	
		// "Hide" the header and footer that we used for the sizing. We need to keep
		// the content of the cell so that the width applied to the header and body
		// both match, but we want to hide it completely. We want to also fix their
		// width to what they currently are
		_fnApplyToChildren( function(nSizer, i) {
			nSizer.innerHTML = '<div class="dataTables_sizing">'+headerContent[i]+'</div>';
			nSizer.childNodes[0].style.height = "0";
			nSizer.childNodes[0].style.overflow = "hidden";
			nSizer.style.width = headerWidths[i];
		}, headerSrcEls );
	
		if ( footer )
		{
			_fnApplyToChildren( function(nSizer, i) {
				nSizer.innerHTML = '<div class="dataTables_sizing">'+footerContent[i]+'</div>';
				nSizer.childNodes[0].style.height = "0";
				nSizer.childNodes[0].style.overflow = "hidden";
				nSizer.style.width = footerWidths[i];
			}, footerSrcEls );
		}
	
		// Sanity check that the table is of a sensible width. If not then we are going to get
		// misalignment - try to prevent this by not allowing the table to shrink below its min width
		if ( Math.round(table.outerWidth()) < Math.round(sanityWidth) )
		{
			// The min width depends upon if we have a vertical scrollbar visible or not */
			correction = ((divBodyEl.scrollHeight > divBodyEl.offsetHeight ||
				divBody.css('overflow-y') == "scroll")) ?
					sanityWidth+barWidth :
					sanityWidth;
	
			// IE6/7 are a law unto themselves...
			if ( ie67 && (divBodyEl.scrollHeight >
				divBodyEl.offsetHeight || divBody.css('overflow-y') == "scroll")
			) {
				tableStyle.width = _fnStringToCss( correction-barWidth );
			}
	
			// And give the user a warning that we've stopped the table getting too small
			if ( scrollX === "" || scrollXInner !== "" ) {
				_fnLog( settings, 1, 'Possible column misalignment', 6 );
			}
		}
		else
		{
			correction = '100%';
		}
	
		// Apply to the container elements
		divBodyStyle.width = _fnStringToCss( correction );
		divHeaderStyle.width = _fnStringToCss( correction );
	
		if ( footer ) {
			settings.nScrollFoot.style.width = _fnStringToCss( correction );
		}
	
	
		/*
		 * 4. Clean up
		 */
		if ( ! scrollY ) {
			/* IE7< puts a vertical scrollbar in place (when it shouldn't be) due to subtracting
			 * the scrollbar height from the visible display, rather than adding it on. We need to
			 * set the height in order to sort this. Don't want to do it in any other browsers.
			 */
			if ( ie67 ) {
				divBodyStyle.height = _fnStringToCss( tableEl.offsetHeight+barWidth );
			}
		}
	
		/* Finally set the width's of the header and footer tables */
		var iOuterWidth = table.outerWidth();
		divHeaderTable[0].style.width = _fnStringToCss( iOuterWidth );
		divHeaderInnerStyle.width = _fnStringToCss( iOuterWidth );
	
		// Figure out if there are scrollbar present - if so then we need a the header and footer to
		// provide a bit more space to allow "overflow" scrolling (i.e. past the scrollbar)
		var bScrolling = table.height() > divBodyEl.clientHeight || divBody.css('overflow-y') == "scroll";
		var padding = 'padding' + (browser.bScrollbarLeft ? 'Left' : 'Right' );
		divHeaderInnerStyle[ padding ] = bScrolling ? barWidth+"px" : "0px";
	
		if ( footer ) {
			divFooterTable[0].style.width = _fnStringToCss( iOuterWidth );
			divFooterInner[0].style.width = _fnStringToCss( iOuterWidth );
			divFooterInner[0].style[padding] = bScrolling ? barWidth+"px" : "0px";
		}
	
		// Correct DOM ordering for colgroup - comes before the thead
		table.children('colgroup').insertBefore( table.children('thead') );
	
		/* Adjust the position of the header in case we loose the y-scrollbar */
		divBody.trigger('scroll');
	
		// If sorting or filtering has occurred, jump the scrolling back to the top
		// only if we aren't holding the position
		if ( (settings.bSorted || settings.bFiltered) && ! settings._drawHold ) {
			divBodyEl.scrollTop = 0;
		}
	}
	
	
	
	/**
	 * Apply a given function to the display child nodes of an element array (typically
	 * TD children of TR rows
	 *  @param {function} fn Method to apply to the objects
	 *  @param array {nodes} an1 List of elements to look through for display children
	 *  @param array {nodes} an2 Another list (identical structure to the first) - optional
	 *  @memberof DataTable#oApi
	 */
	function _fnApplyToChildren( fn, an1, an2 )
	{
		var index=0, i=0, iLen=an1.length;
		var nNode1, nNode2;
	
		while ( i < iLen ) {
			nNode1 = an1[i].firstChild;
			nNode2 = an2 ? an2[i].firstChild : null;
	
			while ( nNode1 ) {
				if ( nNode1.nodeType === 1 ) {
					if ( an2 ) {
						fn( nNode1, nNode2, index );
					}
					else {
						fn( nNode1, index );
					}
	
					index++;
				}
	
				nNode1 = nNode1.nextSibling;
				nNode2 = an2 ? nNode2.nextSibling : null;
			}
	
			i++;
		}
	}
	
	
	
	var __re_html_remove = /<.*?>/g;
	
	
	/**
	 * Calculate the width of columns for the table
	 *  @param {object} oSettings dataTables settings object
	 *  @memberof DataTable#oApi
	 */
	function _fnCalculateColumnWidths ( oSettings )
	{
		var
			table = oSettings.nTable,
			columns = oSettings.aoColumns,
			scroll = oSettings.oScroll,
			scrollY = scroll.sY,
			scrollX = scroll.sX,
			scrollXInner = scroll.sXInner,
			columnCount = columns.length,
			visibleColumns = _fnGetColumns( oSettings, 'bVisible' ),
			headerCells = $('th', oSettings.nTHead),
			tableWidthAttr = table.getAttribute('width'), // from DOM element
			tableContainer = table.parentNode,
			userInputs = false,
			i, column, columnIdx, width, outerWidth,
			browser = oSettings.oBrowser,
			ie67 = browser.bScrollOversize;
	
		var styleWidth = table.style.width;
		if ( styleWidth && styleWidth.indexOf('%') !== -1 ) {
			tableWidthAttr = styleWidth;
		}
	
		/* Convert any user input sizes into pixel sizes */
		for ( i=0 ; i<visibleColumns.length ; i++ ) {
			column = columns[ visibleColumns[i] ];
	
			if ( column.sWidth !== null ) {
				column.sWidth = _fnConvertToWidth( column.sWidthOrig, tableContainer );
	
				userInputs = true;
			}
		}
	
		/* If the number of columns in the DOM equals the number that we have to
		 * process in DataTables, then we can use the offsets that are created by
		 * the web- browser. No custom sizes can be set in order for this to happen,
		 * nor scrolling used
		 */
		if ( ie67 || ! userInputs && ! scrollX && ! scrollY &&
		     columnCount == _fnVisbleColumns( oSettings ) &&
		     columnCount == headerCells.length
		) {
			for ( i=0 ; i<columnCount ; i++ ) {
				var colIdx = _fnVisibleToColumnIndex( oSettings, i );
	
				if ( colIdx !== null ) {
					columns[ colIdx ].sWidth = _fnStringToCss( headerCells.eq(i).width() );
				}
			}
		}
		else
		{
			// Otherwise construct a single row, worst case, table with the widest
			// node in the data, assign any user defined widths, then insert it into
			// the DOM and allow the browser to do all the hard work of calculating
			// table widths
			var tmpTable = $(table).clone() // don't use cloneNode - IE8 will remove events on the main table
				.css( 'visibility', 'hidden' )
				.removeAttr( 'id' );
	
			// Clean up the table body
			tmpTable.find('tbody tr').remove();
			var tr = $('<tr/>').appendTo( tmpTable.find('tbody') );
	
			// Clone the table header and footer - we can't use the header / footer
			// from the cloned table, since if scrolling is active, the table's
			// real header and footer are contained in different table tags
			tmpTable.find('thead, tfoot').remove();
			tmpTable
				.append( $(oSettings.nTHead).clone() )
				.append( $(oSettings.nTFoot).clone() );
	
			// Remove any assigned widths from the footer (from scrolling)
			tmpTable.find('tfoot th, tfoot td').css('width', '');
	
			// Apply custom sizing to the cloned header
			headerCells = _fnGetUniqueThs( oSettings, tmpTable.find('thead')[0] );
	
			for ( i=0 ; i<visibleColumns.length ; i++ ) {
				column = columns[ visibleColumns[i] ];
	
				headerCells[i].style.width = column.sWidthOrig !== null && column.sWidthOrig !== '' ?
					_fnStringToCss( column.sWidthOrig ) :
					'';
	
				// For scrollX we need to force the column width otherwise the
				// browser will collapse it. If this width is smaller than the
				// width the column requires, then it will have no effect
				if ( column.sWidthOrig && scrollX ) {
					$( headerCells[i] ).append( $('<div/>').css( {
						width: column.sWidthOrig,
						margin: 0,
						padding: 0,
						border: 0,
						height: 1
					} ) );
				}
			}
	
			// Find the widest cell for each column and put it into the table
			if ( oSettings.aoData.length ) {
				for ( i=0 ; i<visibleColumns.length ; i++ ) {
					columnIdx = visibleColumns[i];
					column = columns[ columnIdx ];
	
					$( _fnGetWidestNode( oSettings, columnIdx ) )
						.clone( false )
						.append( column.sContentPadding )
						.appendTo( tr );
				}
			}
	
			// Tidy the temporary table - remove name attributes so there aren't
			// duplicated in the dom (radio elements for example)
			$('[name]', tmpTable).removeAttr('name');
	
			// Table has been built, attach to the document so we can work with it.
			// A holding element is used, positioned at the top of the container
			// with minimal height, so it has no effect on if the container scrolls
			// or not. Otherwise it might trigger scrolling when it actually isn't
			// needed
			var holder = $('<div/>').css( scrollX || scrollY ?
					{
						position: 'absolute',
						top: 0,
						left: 0,
						height: 1,
						right: 0,
						overflow: 'hidden'
					} :
					{}
				)
				.append( tmpTable )
				.appendTo( tableContainer );
	
			// When scrolling (X or Y) we want to set the width of the table as 
			// appropriate. However, when not scrolling leave the table width as it
			// is. This results in slightly different, but I think correct behaviour
			if ( scrollX && scrollXInner ) {
				tmpTable.width( scrollXInner );
			}
			else if ( scrollX ) {
				tmpTable.css( 'width', 'auto' );
				tmpTable.removeAttr('width');
	
				// If there is no width attribute or style, then allow the table to
				// collapse
				if ( tmpTable.width() < tableContainer.clientWidth && tableWidthAttr ) {
					tmpTable.width( tableContainer.clientWidth );
				}
			}
			else if ( scrollY ) {
				tmpTable.width( tableContainer.clientWidth );
			}
			else if ( tableWidthAttr ) {
				tmpTable.width( tableWidthAttr );
			}
	
			// Get the width of each column in the constructed table - we need to
			// know the inner width (so it can be assigned to the other table's
			// cells) and the outer width so we can calculate the full width of the
			// table. This is safe since DataTables requires a unique cell for each
			// column, but if ever a header can span multiple columns, this will
			// need to be modified.
			var total = 0;
			for ( i=0 ; i<visibleColumns.length ; i++ ) {
				var cell = $(headerCells[i]);
				var border = cell.outerWidth() - cell.width();
	
				// Use getBounding... where possible (not IE8-) because it can give
				// sub-pixel accuracy, which we then want to round up!
				var bounding = browser.bBounding ?
					Math.ceil( headerCells[i].getBoundingClientRect().width ) :
					cell.outerWidth();
	
				// Total is tracked to remove any sub-pixel errors as the outerWidth
				// of the table might not equal the total given here (IE!).
				total += bounding;
	
				// Width for each column to use
				columns[ visibleColumns[i] ].sWidth = _fnStringToCss( bounding - border );
			}
	
			table.style.width = _fnStringToCss( total );
	
			// Finished with the table - ditch it
			holder.remove();
		}
	
		// If there is a width attr, we want to attach an event listener which
		// allows the table sizing to automatically adjust when the window is
		// resized. Use the width attr rather than CSS, since we can't know if the
		// CSS is a relative value or absolute - DOM read is always px.
		if ( tableWidthAttr ) {
			table.style.width = _fnStringToCss( tableWidthAttr );
		}
	
		if ( (tableWidthAttr || scrollX) && ! oSettings._reszEvt ) {
			var bindResize = function () {
				$(window).on('resize.DT-'+oSettings.sInstance, _fnThrottle( function () {
					_fnAdjustColumnSizing( oSettings );
				} ) );
			};
	
			// IE6/7 will crash if we bind a resize event handler on page load.
			// To be removed in 1.11 which drops IE6/7 support
			if ( ie67 ) {
				setTimeout( bindResize, 1000 );
			}
			else {
				bindResize();
			}
	
			oSettings._reszEvt = true;
		}
	}
	
	
	/**
	 * Throttle the calls to a function. Arguments and context are maintained for
	 * the throttled function
	 *  @param {function} fn Function to be called
	 *  @param {int} [freq=200] call frequency in mS
	 *  @returns {function} wrapped function
	 *  @memberof DataTable#oApi
	 */
	var _fnThrottle = DataTable.util.throttle;
	
	
	/**
	 * Convert a CSS unit width to pixels (e.g. 2em)
	 *  @param {string} width width to be converted
	 *  @param {node} parent parent to get the with for (required for relative widths) - optional
	 *  @returns {int} width in pixels
	 *  @memberof DataTable#oApi
	 */
	function _fnConvertToWidth ( width, parent )
	{
		if ( ! width ) {
			return 0;
		}
	
		var n = $('<div/>')
			.css( 'width', _fnStringToCss( width ) )
			.appendTo( parent || document.body );
	
		var val = n[0].offsetWidth;
		n.remove();
	
		return val;
	}
	
	
	/**
	 * Get the widest node
	 *  @param {object} settings dataTables settings object
	 *  @param {int} colIdx column of interest
	 *  @returns {node} widest table node
	 *  @memberof DataTable#oApi
	 */
	function _fnGetWidestNode( settings, colIdx )
	{
		var idx = _fnGetMaxLenString( settings, colIdx );
		if ( idx < 0 ) {
			return null;
		}
	
		var data = settings.aoData[ idx ];
		return ! data.nTr ? // Might not have been created when deferred rendering
			$('<td/>').html( _fnGetCellData( settings, idx, colIdx, 'display' ) )[0] :
			data.anCells[ colIdx ];
	}
	
	
	/**
	 * Get the maximum strlen for each data column
	 *  @param {object} settings dataTables settings object
	 *  @param {int} colIdx column of interest
	 *  @returns {string} max string length for each column
	 *  @memberof DataTable#oApi
	 */
	function _fnGetMaxLenString( settings, colIdx )
	{
		var s, max=-1, maxIdx = -1;
	
		for ( var i=0, ien=settings.aoData.length ; i<ien ; i++ ) {
			s = _fnGetCellData( settings, i, colIdx, 'display' )+'';
			s = s.replace( __re_html_remove, '' );
			s = s.replace( /&nbsp;/g, ' ' );
	
			if ( s.length > max ) {
				max = s.length;
				maxIdx = i;
			}
		}
	
		return maxIdx;
	}
	
	
	/**
	 * Append a CSS unit (only if required) to a string
	 *  @param {string} value to css-ify
	 *  @returns {string} value with css unit
	 *  @memberof DataTable#oApi
	 */
	function _fnStringToCss( s )
	{
		if ( s === null ) {
			return '0px';
		}
	
		if ( typeof s == 'number' ) {
			return s < 0 ?
				'0px' :
				s+'px';
		}
	
		// Check it has a unit character already
		return s.match(/\d$/) ?
			s+'px' :
			s;
	}
	
	
	
	function _fnSortFlatten ( settings )
	{
		var
			i, iLen, k, kLen,
			aSort = [],
			aiOrig = [],
			aoColumns = settings.aoColumns,
			aDataSort, iCol, sType, srcCol,
			fixed = settings.aaSortingFixed,
			fixedObj = $.isPlainObject( fixed ),
			nestedSort = [],
			add = function ( a ) {
				if ( a.length && ! Array.isArray( a[0] ) ) {
					// 1D array
					nestedSort.push( a );
				}
				else {
					// 2D array
					$.merge( nestedSort, a );
				}
			};
	
		// Build the sort array, with pre-fix and post-fix options if they have been
		// specified
		if ( Array.isArray( fixed ) ) {
			add( fixed );
		}
	
		if ( fixedObj && fixed.pre ) {
			add( fixed.pre );
		}
	
		add( settings.aaSorting );
	
		if (fixedObj && fixed.post ) {
			add( fixed.post );
		}
	
		for ( i=0 ; i<nestedSort.length ; i++ )
		{
			srcCol = nestedSort[i][0];
			aDataSort = aoColumns[ srcCol ].aDataSort;
	
			for ( k=0, kLen=aDataSort.length ; k<kLen ; k++ )
			{
				iCol = aDataSort[k];
				sType = aoColumns[ iCol ].sType || 'string';
	
				if ( nestedSort[i]._idx === undefined ) {
					nestedSort[i]._idx = $.inArray( nestedSort[i][1], aoColumns[iCol].asSorting );
				}
	
				aSort.push( {
					src:       srcCol,
					col:       iCol,
					dir:       nestedSort[i][1],
					index:     nestedSort[i]._idx,
					type:      sType,
					formatter: DataTable.ext.type.order[ sType+"-pre" ]
				} );
			}
		}
	
		return aSort;
	}
	
	/**
	 * Change the order of the table
	 *  @param {object} oSettings dataTables settings object
	 *  @memberof DataTable#oApi
	 *  @todo This really needs split up!
	 */
	function _fnSort ( oSettings )
	{
		var
			i, ien, iLen, j, jLen, k, kLen,
			sDataType, nTh,
			aiOrig = [],
			oExtSort = DataTable.ext.type.order,
			aoData = oSettings.aoData,
			aoColumns = oSettings.aoColumns,
			aDataSort, data, iCol, sType, oSort,
			formatters = 0,
			sortCol,
			displayMaster = oSettings.aiDisplayMaster,
			aSort;
	
		// Resolve any column types that are unknown due to addition or invalidation
		// @todo Can this be moved into a 'data-ready' handler which is called when
		//   data is going to be used in the table?
		_fnColumnTypes( oSettings );
	
		aSort = _fnSortFlatten( oSettings );
	
		for ( i=0, ien=aSort.length ; i<ien ; i++ ) {
			sortCol = aSort[i];
	
			// Track if we can use the fast sort algorithm
			if ( sortCol.formatter ) {
				formatters++;
			}
	
			// Load the data needed for the sort, for each cell
			_fnSortData( oSettings, sortCol.col );
		}
	
		/* No sorting required if server-side or no sorting array */
		if ( _fnDataSource( oSettings ) != 'ssp' && aSort.length !== 0 )
		{
			// Create a value - key array of the current row positions such that we can use their
			// current position during the sort, if values match, in order to perform stable sorting
			for ( i=0, iLen=displayMaster.length ; i<iLen ; i++ ) {
				aiOrig[ displayMaster[i] ] = i;
			}
	
			/* Do the sort - here we want multi-column sorting based on a given data source (column)
			 * and sorting function (from oSort) in a certain direction. It's reasonably complex to
			 * follow on it's own, but this is what we want (example two column sorting):
			 *  fnLocalSorting = function(a,b){
			 *    var iTest;
			 *    iTest = oSort['string-asc']('data11', 'data12');
			 *      if (iTest !== 0)
			 *        return iTest;
			 *    iTest = oSort['numeric-desc']('data21', 'data22');
			 *    if (iTest !== 0)
			 *      return iTest;
			 *    return oSort['numeric-asc']( aiOrig[a], aiOrig[b] );
			 *  }
			 * Basically we have a test for each sorting column, if the data in that column is equal,
			 * test the next column. If all columns match, then we use a numeric sort on the row
			 * positions in the original data array to provide a stable sort.
			 *
			 * Note - I know it seems excessive to have two sorting methods, but the first is around
			 * 15% faster, so the second is only maintained for backwards compatibility with sorting
			 * methods which do not have a pre-sort formatting function.
			 */
			if ( formatters === aSort.length ) {
				// All sort types have formatting functions
				displayMaster.sort( function ( a, b ) {
					var
						x, y, k, test, sort,
						len=aSort.length,
						dataA = aoData[a]._aSortData,
						dataB = aoData[b]._aSortData;
	
					for ( k=0 ; k<len ; k++ ) {
						sort = aSort[k];
	
						x = dataA[ sort.col ];
						y = dataB[ sort.col ];
	
						test = x<y ? -1 : x>y ? 1 : 0;
						if ( test !== 0 ) {
							return sort.dir === 'asc' ? test : -test;
						}
					}
	
					x = aiOrig[a];
					y = aiOrig[b];
					return x<y ? -1 : x>y ? 1 : 0;
				} );
			}
			else {
				// Depreciated - remove in 1.11 (providing a plug-in option)
				// Not all sort types have formatting methods, so we have to call their sorting
				// methods.
				displayMaster.sort( function ( a, b ) {
					var
						x, y, k, l, test, sort, fn,
						len=aSort.length,
						dataA = aoData[a]._aSortData,
						dataB = aoData[b]._aSortData;
	
					for ( k=0 ; k<len ; k++ ) {
						sort = aSort[k];
	
						x = dataA[ sort.col ];
						y = dataB[ sort.col ];
	
						fn = oExtSort[ sort.type+"-"+sort.dir ] || oExtSort[ "string-"+sort.dir ];
						test = fn( x, y );
						if ( test !== 0 ) {
							return test;
						}
					}
	
					x = aiOrig[a];
					y = aiOrig[b];
					return x<y ? -1 : x>y ? 1 : 0;
				} );
			}
		}
	
		/* Tell the draw function that we have sorted the data */
		oSettings.bSorted = true;
	}
	
	
	function _fnSortAria ( settings )
	{
		var label;
		var nextSort;
		var columns = settings.aoColumns;
		var aSort = _fnSortFlatten( settings );
		var oAria = settings.oLanguage.oAria;
	
		// ARIA attributes - need to loop all columns, to update all (removing old
		// attributes as needed)
		for ( var i=0, iLen=columns.length ; i<iLen ; i++ )
		{
			var col = columns[i];
			var asSorting = col.asSorting;
			var sTitle = col.ariaTitle || col.sTitle.replace( /<.*?>/g, "" );
			var th = col.nTh;
	
			// IE7 is throwing an error when setting these properties with jQuery's
			// attr() and removeAttr() methods...
			th.removeAttribute('aria-sort');
	
			/* In ARIA only the first sorting column can be marked as sorting - no multi-sort option */
			if ( col.bSortable ) {
				if ( aSort.length > 0 && aSort[0].col == i ) {
					th.setAttribute('aria-sort', aSort[0].dir=="asc" ? "ascending" : "descending" );
					nextSort = asSorting[ aSort[0].index+1 ] || asSorting[0];
				}
				else {
					nextSort = asSorting[0];
				}
	
				label = sTitle + ( nextSort === "asc" ?
					oAria.sSortAscending :
					oAria.sSortDescending
				);
			}
			else {
				label = sTitle;
			}
	
			th.setAttribute('aria-label', label);
		}
	}
	
	
	/**
	 * Function to run on user sort request
	 *  @param {object} settings dataTables settings object
	 *  @param {node} attachTo node to attach the handler to
	 *  @param {int} colIdx column sorting index
	 *  @param {boolean} [append=false] Append the requested sort to the existing
	 *    sort if true (i.e. multi-column sort)
	 *  @param {function} [callback] callback function
	 *  @memberof DataTable#oApi
	 */
	function _fnSortListener ( settings, colIdx, append, callback )
	{
		var col = settings.aoColumns[ colIdx ];
		var sorting = settings.aaSorting;
		var asSorting = col.asSorting;
		var nextSortIdx;
		var next = function ( a, overflow ) {
			var idx = a._idx;
			if ( idx === undefined ) {
				idx = $.inArray( a[1], asSorting );
			}
	
			return idx+1 < asSorting.length ?
				idx+1 :
				overflow ?
					null :
					0;
		};
	
		// Convert to 2D array if needed
		if ( typeof sorting[0] === 'number' ) {
			sorting = settings.aaSorting = [ sorting ];
		}
	
		// If appending the sort then we are multi-column sorting
		if ( append && settings.oFeatures.bSortMulti ) {
			// Are we already doing some kind of sort on this column?
			var sortIdx = $.inArray( colIdx, _pluck(sorting, '0') );
	
			if ( sortIdx !== -1 ) {
				// Yes, modify the sort
				nextSortIdx = next( sorting[sortIdx], true );
	
				if ( nextSortIdx === null && sorting.length === 1 ) {
					nextSortIdx = 0; // can't remove sorting completely
				}
	
				if ( nextSortIdx === null ) {
					sorting.splice( sortIdx, 1 );
				}
				else {
					sorting[sortIdx][1] = asSorting[ nextSortIdx ];
					sorting[sortIdx]._idx = nextSortIdx;
				}
			}
			else {
				// No sort on this column yet
				sorting.push( [ colIdx, asSorting[0], 0 ] );
				sorting[sorting.length-1]._idx = 0;
			}
		}
		else if ( sorting.length && sorting[0][0] == colIdx ) {
			// Single column - already sorting on this column, modify the sort
			nextSortIdx = next( sorting[0] );
	
			sorting.length = 1;
			sorting[0][1] = asSorting[ nextSortIdx ];
			sorting[0]._idx = nextSortIdx;
		}
		else {
			// Single column - sort only on this column
			sorting.length = 0;
			sorting.push( [ colIdx, asSorting[0] ] );
			sorting[0]._idx = 0;
		}
	
		// Run the sort by calling a full redraw
		_fnReDraw( settings );
	
		// callback used for async user interaction
		if ( typeof callback == 'function' ) {
			callback( settings );
		}
	}
	
	
	/**
	 * Attach a sort handler (click) to a node
	 *  @param {object} settings dataTables settings object
	 *  @param {node} attachTo node to attach the handler to
	 *  @param {int} colIdx column sorting index
	 *  @param {function} [callback] callback function
	 *  @memberof DataTable#oApi
	 */
	function _fnSortAttachListener ( settings, attachTo, colIdx, callback )
	{
		var col = settings.aoColumns[ colIdx ];
	
		_fnBindAction( attachTo, {}, function (e) {
			/* If the column is not sortable - don't to anything */
			if ( col.bSortable === false ) {
				return;
			}
	
			// If processing is enabled use a timeout to allow the processing
			// display to be shown - otherwise to it synchronously
			if ( settings.oFeatures.bProcessing ) {
				_fnProcessingDisplay( settings, true );
	
				setTimeout( function() {
					_fnSortListener( settings, colIdx, e.shiftKey, callback );
	
					// In server-side processing, the draw callback will remove the
					// processing display
					if ( _fnDataSource( settings ) !== 'ssp' ) {
						_fnProcessingDisplay( settings, false );
					}
				}, 0 );
			}
			else {
				_fnSortListener( settings, colIdx, e.shiftKey, callback );
			}
		} );
	}
	
	
	/**
	 * Set the sorting classes on table's body, Note: it is safe to call this function
	 * when bSort and bSortClasses are false
	 *  @param {object} oSettings dataTables settings object
	 *  @memberof DataTable#oApi
	 */
	function _fnSortingClasses( settings )
	{
		var oldSort = settings.aLastSort;
		var sortClass = settings.oClasses.sSortColumn;
		var sort = _fnSortFlatten( settings );
		var features = settings.oFeatures;
		var i, ien, colIdx;
	
		if ( features.bSort && features.bSortClasses ) {
			// Remove old sorting classes
			for ( i=0, ien=oldSort.length ; i<ien ; i++ ) {
				colIdx = oldSort[i].src;
	
				// Remove column sorting
				$( _pluck( settings.aoData, 'anCells', colIdx ) )
					.removeClass( sortClass + (i<2 ? i+1 : 3) );
			}
	
			// Add new column sorting
			for ( i=0, ien=sort.length ; i<ien ; i++ ) {
				colIdx = sort[i].src;
	
				$( _pluck( settings.aoData, 'anCells', colIdx ) )
					.addClass( sortClass + (i<2 ? i+1 : 3) );
			}
		}
	
		settings.aLastSort = sort;
	}
	
	
	// Get the data to sort a column, be it from cache, fresh (populating the
	// cache), or from a sort formatter
	function _fnSortData( settings, idx )
	{
		// Custom sorting function - provided by the sort data type
		var column = settings.aoColumns[ idx ];
		var customSort = DataTable.ext.order[ column.sSortDataType ];
		var customData;
	
		if ( customSort ) {
			customData = customSort.call( settings.oInstance, settings, idx,
				_fnColumnIndexToVisible( settings, idx )
			);
		}
	
		// Use / populate cache
		var row, cellData;
		var formatter = DataTable.ext.type.order[ column.sType+"-pre" ];
	
		for ( var i=0, ien=settings.aoData.length ; i<ien ; i++ ) {
			row = settings.aoData[i];
	
			if ( ! row._aSortData ) {
				row._aSortData = [];
			}
	
			if ( ! row._aSortData[idx] || customSort ) {
				cellData = customSort ?
					customData[i] : // If there was a custom sort function, use data from there
					_fnGetCellData( settings, i, idx, 'sort' );
	
				row._aSortData[ idx ] = formatter ?
					formatter( cellData ) :
					cellData;
			}
		}
	}
	
	
	
	/**
	 * Save the state of a table
	 *  @param {object} oSettings dataTables settings object
	 *  @memberof DataTable#oApi
	 */
	function _fnSaveState ( settings )
	{
		if (settings._bLoadingState) {
			return;
		}
	
		/* Store the interesting variables */
		var state = {
			time:    +new Date(),
			start:   settings._iDisplayStart,
			length:  settings._iDisplayLength,
			order:   $.extend( true, [], settings.aaSorting ),
			search:  _fnSearchToCamel( settings.oPreviousSearch ),
			columns: $.map( settings.aoColumns, function ( col, i ) {
				return {
					visible: col.bVisible,
					search: _fnSearchToCamel( settings.aoPreSearchCols[i] )
				};
			} )
		};
	
		settings.oSavedState = state;
		_fnCallbackFire( settings, "aoStateSaveParams", 'stateSaveParams', [settings, state] );
		
		if ( settings.oFeatures.bStateSave && !settings.bDestroying )
		{
			settings.fnStateSaveCallback.call( settings.oInstance, settings, state );
		}	
	}
	
	
	/**
	 * Attempt to load a saved table state
	 *  @param {object} oSettings dataTables settings object
	 *  @param {object} oInit DataTables init object so we can override settings
	 *  @param {function} callback Callback to execute when the state has been loaded
	 *  @memberof DataTable#oApi
	 */
	function _fnLoadState ( settings, oInit, callback )
	{
		if ( ! settings.oFeatures.bStateSave ) {
			callback();
			return;
		}
	
		var loaded = function(state) {
			_fnImplementState(settings, state, callback);
		}
	
		var state = settings.fnStateLoadCallback.call( settings.oInstance, settings, loaded );
	
		if ( state !== undefined ) {
			_fnImplementState( settings, state, callback );
		}
		// otherwise, wait for the loaded callback to be executed
	
		return true;
	}
	
	function _fnImplementState ( settings, s, callback) {
		var i, ien;
		var columns = settings.aoColumns;
		settings._bLoadingState = true;
	
		// When StateRestore was introduced the state could now be implemented at any time
		// Not just initialisation. To do this an api instance is required in some places
		var api = settings._bInitComplete ? new DataTable.Api(settings) : null;
	
		if ( ! s || ! s.time ) {
			settings._bLoadingState = false;
			callback();
			return;
		}
	
		// Allow custom and plug-in manipulation functions to alter the saved data set and
		// cancelling of loading by returning false
		var abStateLoad = _fnCallbackFire( settings, 'aoStateLoadParams', 'stateLoadParams', [settings, s] );
		if ( $.inArray( false, abStateLoad ) !== -1 ) {
			settings._bLoadingState = false;
			callback();
			return;
		}
	
		// Reject old data
		var duration = settings.iStateDuration;
		if ( duration > 0 && s.time < +new Date() - (duration*1000) ) {
			settings._bLoadingState = false;
			callback();
			return;
		}
	
		// Number of columns have changed - all bets are off, no restore of settings
		if ( s.columns && columns.length !== s.columns.length ) {
			settings._bLoadingState = false;
			callback();
			return;
		}
	
		// Store the saved state so it might be accessed at any time
		settings.oLoadedState = $.extend( true, {}, s );
	
		// Page Length
		if ( s.length !== undefined ) {
			// If already initialised just set the value directly so that the select element is also updated
			if (api) {
				api.page.len(s.length)
			}
			else {
				settings._iDisplayLength   = s.length;
			}
		}
	
		// Restore key features - todo - for 1.11 this needs to be done by
		// subscribed events
		if ( s.start !== undefined ) {
			if(api === null) {
				settings._iDisplayStart    = s.start;
				settings.iInitDisplayStart = s.start;
			}
			else {
				_fnPageChange(settings, s.start/settings._iDisplayLength);
			}
		}
	
		// Order
		if ( s.order !== undefined ) {
			settings.aaSorting = [];
			$.each( s.order, function ( i, col ) {
				settings.aaSorting.push( col[0] >= columns.length ?
					[ 0, col[1] ] :
					col
				);
			} );
		}
	
		// Search
		if ( s.search !== undefined ) {
			$.extend( settings.oPreviousSearch, _fnSearchToHung( s.search ) );
		}
	
		// Columns
		if ( s.columns ) {
			for ( i=0, ien=s.columns.length ; i<ien ; i++ ) {
				var col = s.columns[i];
	
				// Visibility
				if ( col.visible !== undefined ) {
					// If the api is defined, the table has been initialised so we need to use it rather than internal settings
					if (api) {
						// Don't redraw the columns on every iteration of this loop, we will do this at the end instead
						api.column(i).visible(col.visible, false);
					}
					else {
						columns[i].bVisible = col.visible;
					}
				}
	
				// Search
				if ( col.search !== undefined ) {
					$.extend( settings.aoPreSearchCols[i], _fnSearchToHung( col.search ) );
				}
			}
			
			// If the api is defined then we need to adjust the columns once the visibility has been changed
			if (api) {
				api.columns.adjust();
			}
		}
	
		settings._bLoadingState = false;
		_fnCallbackFire( settings, 'aoStateLoaded', 'stateLoaded', [settings, s] );
		callback();
	};
	
	
	/**
	 * Return the settings object for a particular table
	 *  @param {node} table table we are using as a dataTable
	 *  @returns {object} Settings object - or null if not found
	 *  @memberof DataTable#oApi
	 */
	function _fnSettingsFromNode ( table )
	{
		var settings = DataTable.settings;
		var idx = $.inArray( table, _pluck( settings, 'nTable' ) );
	
		return idx !== -1 ?
			settings[ idx ] :
			null;
	}
	
	
	/**
	 * Log an error message
	 *  @param {object} settings dataTables settings object
	 *  @param {int} level log error messages, or display them to the user
	 *  @param {string} msg error message
	 *  @param {int} tn Technical note id to get more information about the error.
	 *  @memberof DataTable#oApi
	 */
	function _fnLog( settings, level, msg, tn )
	{
		msg = 'DataTables warning: '+
			(settings ? 'table id='+settings.sTableId+' - ' : '')+msg;
	
		if ( tn ) {
			msg += '. For more information about this error, please see '+
			'http://datatables.net/tn/'+tn;
		}
	
		if ( ! level  ) {
			// Backwards compatibility pre 1.10
			var ext = DataTable.ext;
			var type = ext.sErrMode || ext.errMode;
	
			if ( settings ) {
				_fnCallbackFire( settings, null, 'error', [ settings, tn, msg ] );
			}
	
			if ( type == 'alert' ) {
				alert( msg );
			}
			else if ( type == 'throw' ) {
				throw new Error(msg);
			}
			else if ( typeof type == 'function' ) {
				type( settings, tn, msg );
			}
		}
		else if ( window.console && console.log ) {
			console.log( msg );
		}
	}
	
	
	/**
	 * See if a property is defined on one object, if so assign it to the other object
	 *  @param {object} ret target object
	 *  @param {object} src source object
	 *  @param {string} name property
	 *  @param {string} [mappedName] name to map too - optional, name used if not given
	 *  @memberof DataTable#oApi
	 */
	function _fnMap( ret, src, name, mappedName )
	{
		if ( Array.isArray( name ) ) {
			$.each( name, function (i, val) {
				if ( Array.isArray( val ) ) {
					_fnMap( ret, src, val[0], val[1] );
				}
				else {
					_fnMap( ret, src, val );
				}
			} );
	
			return;
		}
	
		if ( mappedName === undefined ) {
			mappedName = name;
		}
	
		if ( src[name] !== undefined ) {
			ret[mappedName] = src[name];
		}
	}
	
	
	/**
	 * Extend objects - very similar to jQuery.extend, but deep copy objects, and
	 * shallow copy arrays. The reason we need to do this, is that we don't want to
	 * deep copy array init values (such as aaSorting) since the dev wouldn't be
	 * able to override them, but we do want to deep copy arrays.
	 *  @param {object} out Object to extend
	 *  @param {object} extender Object from which the properties will be applied to
	 *      out
	 *  @param {boolean} breakRefs If true, then arrays will be sliced to take an
	 *      independent copy with the exception of the `data` or `aaData` parameters
	 *      if they are present. This is so you can pass in a collection to
	 *      DataTables and have that used as your data source without breaking the
	 *      references
	 *  @returns {object} out Reference, just for convenience - out === the return.
	 *  @memberof DataTable#oApi
	 *  @todo This doesn't take account of arrays inside the deep copied objects.
	 */
	function _fnExtend( out, extender, breakRefs )
	{
		var val;
	
		for ( var prop in extender ) {
			if ( extender.hasOwnProperty(prop) ) {
				val = extender[prop];
	
				if ( $.isPlainObject( val ) ) {
					if ( ! $.isPlainObject( out[prop] ) ) {
						out[prop] = {};
					}
					$.extend( true, out[prop], val );
				}
				else if ( breakRefs && prop !== 'data' && prop !== 'aaData' && Array.isArray(val) ) {
					out[prop] = val.slice();
				}
				else {
					out[prop] = val;
				}
			}
		}
	
		return out;
	}
	
	
	/**
	 * Bind an event handers to allow a click or return key to activate the callback.
	 * This is good for accessibility since a return on the keyboard will have the
	 * same effect as a click, if the element has focus.
	 *  @param {element} n Element to bind the action to
	 *  @param {object} oData Data object to pass to the triggered function
	 *  @param {function} fn Callback function for when the event is triggered
	 *  @memberof DataTable#oApi
	 */
	function _fnBindAction( n, oData, fn )
	{
		$(n)
			.on( 'click.DT', oData, function (e) {
					$(n).trigger('blur'); // Remove focus outline for mouse users
					fn(e);
				} )
			.on( 'keypress.DT', oData, function (e){
					if ( e.which === 13 ) {
						e.preventDefault();
						fn(e);
					}
				} )
			.on( 'selectstart.DT', function () {
					/* Take the brutal approach to cancelling text selection */
					return false;
				} );
	}
	
	
	/**
	 * Register a callback function. Easily allows a callback function to be added to
	 * an array store of callback functions that can then all be called together.
	 *  @param {object} oSettings dataTables settings object
	 *  @param {string} sStore Name of the array storage for the callbacks in oSettings
	 *  @param {function} fn Function to be called back
	 *  @param {string} sName Identifying name for the callback (i.e. a label)
	 *  @memberof DataTable#oApi
	 */
	function _fnCallbackReg( oSettings, sStore, fn, sName )
	{
		if ( fn )
		{
			oSettings[sStore].push( {
				"fn": fn,
				"sName": sName
			} );
		}
	}
	
	
	/**
	 * Fire callback functions and trigger events. Note that the loop over the
	 * callback array store is done backwards! Further note that you do not want to
	 * fire off triggers in time sensitive applications (for example cell creation)
	 * as its slow.
	 *  @param {object} settings dataTables settings object
	 *  @param {string} callbackArr Name of the array storage for the callbacks in
	 *      oSettings
	 *  @param {string} eventName Name of the jQuery custom event to trigger. If
	 *      null no trigger is fired
	 *  @param {array} args Array of arguments to pass to the callback function /
	 *      trigger
	 *  @memberof DataTable#oApi
	 */
	function _fnCallbackFire( settings, callbackArr, eventName, args )
	{
		var ret = [];
	
		if ( callbackArr ) {
			ret = $.map( settings[callbackArr].slice().reverse(), function (val, i) {
				return val.fn.apply( settings.oInstance, args );
			} );
		}
	
		if ( eventName !== null ) {
			var e = $.Event( eventName+'.dt' );
			var table = $(settings.nTable);
	
			table.trigger( e, args );
	
			// If not yet attached to the document, trigger the event
			// on the body directly to sort of simulate the bubble
			if (table.parents('body').length === 0) {
				$('body').trigger( e, args );
			}
	
			ret.push( e.result );
		}
	
		return ret;
	}
	
	
	function _fnLengthOverflow ( settings )
	{
		var
			start = settings._iDisplayStart,
			end = settings.fnDisplayEnd(),
			len = settings._iDisplayLength;
	
		/* If we have space to show extra rows (backing up from the end point - then do so */
		if ( start >= end )
		{
			start = end - len;
		}
	
		// Keep the start record on the current page
		start -= (start % len);
	
		if ( len === -1 || start < 0 )
		{
			start = 0;
		}
	
		settings._iDisplayStart = start;
	}
	
	
	function _fnRenderer( settings, type )
	{
		var renderer = settings.renderer;
		var host = DataTable.ext.renderer[type];
	
		if ( $.isPlainObject( renderer ) && renderer[type] ) {
			// Specific renderer for this type. If available use it, otherwise use
			// the default.
			return host[renderer[type]] || host._;
		}
		else if ( typeof renderer === 'string' ) {
			// Common renderer - if there is one available for this type use it,
			// otherwise use the default
			return host[renderer] || host._;
		}
	
		// Use the default
		return host._;
	}
	
	
	/**
	 * Detect the data source being used for the table. Used to simplify the code
	 * a little (ajax) and to make it compress a little smaller.
	 *
	 *  @param {object} settings dataTables settings object
	 *  @returns {string} Data source
	 *  @memberof DataTable#oApi
	 */
	function _fnDataSource ( settings )
	{
		if ( settings.oFeatures.bServerSide ) {
			return 'ssp';
		}
		else if ( settings.ajax || settings.sAjaxSource ) {
			return 'ajax';
		}
		return 'dom';
	}
	
	
	
	
	/**
	 * Computed structure of the DataTables API, defined by the options passed to
	 * `DataTable.Api.register()` when building the API.
	 *
	 * The structure is built in order to speed creation and extension of the Api
	 * objects since the extensions are effectively pre-parsed.
	 *
	 * The array is an array of objects with the following structure, where this
	 * base array represents the Api prototype base:
	 *
	 *     [
	 *       {
	 *         name:      'data'                -- string   - Property name
	 *         val:       function () {},       -- function - Api method (or undefined if just an object
	 *         methodExt: [ ... ],              -- array    - Array of Api object definitions to extend the method result
	 *         propExt:   [ ... ]               -- array    - Array of Api object definitions to extend the property
	 *       },
	 *       {
	 *         name:     'row'
	 *         val:       {},
	 *         methodExt: [ ... ],
	 *         propExt:   [
	 *           {
	 *             name:      'data'
	 *             val:       function () {},
	 *             methodExt: [ ... ],
	 *             propExt:   [ ... ]
	 *           },
	 *           ...
	 *         ]
	 *       }
	 *     ]
	 *
	 * @type {Array}
	 * @ignore
	 */
	var __apiStruct = [];
	
	
	/**
	 * `Array.prototype` reference.
	 *
	 * @type object
	 * @ignore
	 */
	var __arrayProto = Array.prototype;
	
	
	/**
	 * Abstraction for `context` parameter of the `Api` constructor to allow it to
	 * take several different forms for ease of use.
	 *
	 * Each of the input parameter types will be converted to a DataTables settings
	 * object where possible.
	 *
	 * @param  {string|node|jQuery|object} mixed DataTable identifier. Can be one
	 *   of:
	 *
	 *   * `string` - jQuery selector. Any DataTables' matching the given selector
	 *     with be found and used.
	 *   * `node` - `TABLE` node which has already been formed into a DataTable.
	 *   * `jQuery` - A jQuery object of `TABLE` nodes.
	 *   * `object` - DataTables settings object
	 *   * `DataTables.Api` - API instance
	 * @return {array|null} Matching DataTables settings objects. `null` or
	 *   `undefined` is returned if no matching DataTable is found.
	 * @ignore
	 */
	var _toSettings = function ( mixed )
	{
		var idx, jq;
		var settings = DataTable.settings;
		var tables = $.map( settings, function (el, i) {
			return el.nTable;
		} );
	
		if ( ! mixed ) {
			return [];
		}
		else if ( mixed.nTable && mixed.oApi ) {
			// DataTables settings object
			return [ mixed ];
		}
		else if ( mixed.nodeName && mixed.nodeName.toLowerCase() === 'table' ) {
			// Table node
			idx = $.inArray( mixed, tables );
			return idx !== -1 ? [ settings[idx] ] : null;
		}
		else if ( mixed && typeof mixed.settings === 'function' ) {
			return mixed.settings().toArray();
		}
		else if ( typeof mixed === 'string' ) {
			// jQuery selector
			jq = $(mixed);
		}
		else if ( mixed instanceof $ ) {
			// jQuery object (also DataTables instance)
			jq = mixed;
		}
	
		if ( jq ) {
			return jq.map( function(i) {
				idx = $.inArray( this, tables );
				return idx !== -1 ? settings[idx] : null;
			} ).toArray();
		}
	};
	
	
	/**
	 * DataTables API class - used to control and interface with  one or more
	 * DataTables enhanced tables.
	 *
	 * The API class is heavily based on jQuery, presenting a chainable interface
	 * that you can use to interact with tables. Each instance of the API class has
	 * a "context" - i.e. the tables that it will operate on. This could be a single
	 * table, all tables on a page or a sub-set thereof.
	 *
	 * Additionally the API is designed to allow you to easily work with the data in
	 * the tables, retrieving and manipulating it as required. This is done by
	 * presenting the API class as an array like interface. The contents of the
	 * array depend upon the actions requested by each method (for example
	 * `rows().nodes()` will return an array of nodes, while `rows().data()` will
	 * return an array of objects or arrays depending upon your table's
	 * configuration). The API object has a number of array like methods (`push`,
	 * `pop`, `reverse` etc) as well as additional helper methods (`each`, `pluck`,
	 * `unique` etc) to assist your working with the data held in a table.
	 *
	 * Most methods (those which return an Api instance) are chainable, which means
	 * the return from a method call also has all of the methods available that the
	 * top level object had. For example, these two calls are equivalent:
	 *
	 *     // Not chained
	 *     api.row.add( {...} );
	 *     api.draw();
	 *
	 *     // Chained
	 *     api.row.add( {...} ).draw();
	 *
	 * @class DataTable.Api
	 * @param {array|object|string|jQuery} context DataTable identifier. This is
	 *   used to define which DataTables enhanced tables this API will operate on.
	 *   Can be one of:
	 *
	 *   * `string` - jQuery selector. Any DataTables' matching the given selector
	 *     with be found and used.
	 *   * `node` - `TABLE` node which has already been formed into a DataTable.
	 *   * `jQuery` - A jQuery object of `TABLE` nodes.
	 *   * `object` - DataTables settings object
	 * @param {array} [data] Data to initialise the Api instance with.
	 *
	 * @example
	 *   // Direct initialisation during DataTables construction
	 *   var api = $('#example').DataTable();
	 *
	 * @example
	 *   // Initialisation using a DataTables jQuery object
	 *   var api = $('#example').dataTable().api();
	 *
	 * @example
	 *   // Initialisation as a constructor
	 *   var api = new $.fn.DataTable.Api( 'table.dataTable' );
	 */
	_Api = function ( context, data )
	{
		if ( ! (this instanceof _Api) ) {
			return new _Api( context, data );
		}
	
		var settings = [];
		var ctxSettings = function ( o ) {
			var a = _toSettings( o );
			if ( a ) {
				settings.push.apply( settings, a );
			}
		};
	
		if ( Array.isArray( context ) ) {
			for ( var i=0, ien=context.length ; i<ien ; i++ ) {
				ctxSettings( context[i] );
			}
		}
		else {
			ctxSettings( context );
		}
	
		// Remove duplicates
		this.context = _unique( settings );
	
		// Initial data
		if ( data ) {
			$.merge( this, data );
		}
	
		// selector
		this.selector = {
			rows: null,
			cols: null,
			opts: null
		};
	
		_Api.extend( this, this, __apiStruct );
	};
	
	DataTable.Api = _Api;
	
	// Don't destroy the existing prototype, just extend it. Required for jQuery 2's
	// isPlainObject.
	$.extend( _Api.prototype, {
		any: function ()
		{
			return this.count() !== 0;
		},
	
	
		concat:  __arrayProto.concat,
	
	
		context: [], // array of table settings objects
	
	
		count: function ()
		{
			return this.flatten().length;
		},
	
	
		each: function ( fn )
		{
			for ( var i=0, ien=this.length ; i<ien; i++ ) {
				fn.call( this, this[i], i, this );
			}
	
			return this;
		},
	
	
		eq: function ( idx )
		{
			var ctx = this.context;
	
			return ctx.length > idx ?
				new _Api( ctx[idx], this[idx] ) :
				null;
		},
	
	
		filter: function ( fn )
		{
			var a = [];
	
			if ( __arrayProto.filter ) {
				a = __arrayProto.filter.call( this, fn, this );
			}
			else {
				// Compatibility for browsers without EMCA-252-5 (JS 1.6)
				for ( var i=0, ien=this.length ; i<ien ; i++ ) {
					if ( fn.call( this, this[i], i, this ) ) {
						a.push( this[i] );
					}
				}
			}
	
			return new _Api( this.context, a );
		},
	
	
		flatten: function ()
		{
			var a = [];
			return new _Api( this.context, a.concat.apply( a, this.toArray() ) );
		},
	
	
		join:    __arrayProto.join,
	
	
		indexOf: __arrayProto.indexOf || function (obj, start)
		{
			for ( var i=(start || 0), ien=this.length ; i<ien ; i++ ) {
				if ( this[i] === obj ) {
					return i;
				}
			}
			return -1;
		},
	
		iterator: function ( flatten, type, fn, alwaysNew ) {
			var
				a = [], ret,
				i, ien, j, jen,
				context = this.context,
				rows, items, item,
				selector = this.selector;
	
			// Argument shifting
			if ( typeof flatten === 'string' ) {
				alwaysNew = fn;
				fn = type;
				type = flatten;
				flatten = false;
			}
	
			for ( i=0, ien=context.length ; i<ien ; i++ ) {
				var apiInst = new _Api( context[i] );
	
				if ( type === 'table' ) {
					ret = fn.call( apiInst, context[i], i );
	
					if ( ret !== undefined ) {
						a.push( ret );
					}
				}
				else if ( type === 'columns' || type === 'rows' ) {
					// this has same length as context - one entry for each table
					ret = fn.call( apiInst, context[i], this[i], i );
	
					if ( ret !== undefined ) {
						a.push( ret );
					}
				}
				else if ( type === 'column' || type === 'column-rows' || type === 'row' || type === 'cell' ) {
					// columns and rows share the same structure.
					// 'this' is an array of column indexes for each context
					items = this[i];
	
					if ( type === 'column-rows' ) {
						rows = _selector_row_indexes( context[i], selector.opts );
					}
	
					for ( j=0, jen=items.length ; j<jen ; j++ ) {
						item = items[j];
	
						if ( type === 'cell' ) {
							ret = fn.call( apiInst, context[i], item.row, item.column, i, j );
						}
						else {
							ret = fn.call( apiInst, context[i], item, i, j, rows );
						}
	
						if ( ret !== undefined ) {
							a.push( ret );
						}
					}
				}
			}
	
			if ( a.length || alwaysNew ) {
				var api = new _Api( context, flatten ? a.concat.apply( [], a ) : a );
				var apiSelector = api.selector;
				apiSelector.rows = selector.rows;
				apiSelector.cols = selector.cols;
				apiSelector.opts = selector.opts;
				return api;
			}
			return this;
		},
	
	
		lastIndexOf: __arrayProto.lastIndexOf || function (obj, start)
		{
			// Bit cheeky...
			return this.indexOf.apply( this.toArray.reverse(), arguments );
		},
	
	
		length:  0,
	
	
		map: function ( fn )
		{
			var a = [];
	
			if ( __arrayProto.map ) {
				a = __arrayProto.map.call( this, fn, this );
			}
			else {
				// Compatibility for browsers without EMCA-252-5 (JS 1.6)
				for ( var i=0, ien=this.length ; i<ien ; i++ ) {
					a.push( fn.call( this, this[i], i ) );
				}
			}
	
			return new _Api( this.context, a );
		},
	
	
		pluck: function ( prop )
		{
			var fn = DataTable.util.get(prop);
	
			return this.map( function ( el ) {
				return fn(el);
			} );
		},
	
		pop:     __arrayProto.pop,
	
	
		push:    __arrayProto.push,
	
	
		// Does not return an API instance
		reduce: __arrayProto.reduce || function ( fn, init )
		{
			return _fnReduce( this, fn, init, 0, this.length, 1 );
		},
	
	
		reduceRight: __arrayProto.reduceRight || function ( fn, init )
		{
			return _fnReduce( this, fn, init, this.length-1, -1, -1 );
		},
	
	
		reverse: __arrayProto.reverse,
	
	
		// Object with rows, columns and opts
		selector: null,
	
	
		shift:   __arrayProto.shift,
	
	
		slice: function () {
			return new _Api( this.context, this );
		},
	
	
		sort:    __arrayProto.sort, // ? name - order?
	
	
		splice:  __arrayProto.splice,
	
	
		toArray: function ()
		{
			return __arrayProto.slice.call( this );
		},
	
	
		to$: function ()
		{
			return $( this );
		},
	
	
		toJQuery: function ()
		{
			return $( this );
		},
	
	
		unique: function ()
		{
			return new _Api( this.context, _unique(this) );
		},
	
	
		unshift: __arrayProto.unshift
	} );
	
	
	_Api.extend = function ( scope, obj, ext )
	{
		// Only extend API instances and static properties of the API
		if ( ! ext.length || ! obj || ( ! (obj instanceof _Api) && ! obj.__dt_wrapper ) ) {
			return;
		}
	
		var
			i, ien,
			struct,
			methodScoping = function ( scope, fn, struc ) {
				return function () {
					var ret = fn.apply( scope, arguments );
	
					// Method extension
					_Api.extend( ret, ret, struc.methodExt );
					return ret;
				};
			};
	
		for ( i=0, ien=ext.length ; i<ien ; i++ ) {
			struct = ext[i];
	
			// Value
			obj[ struct.name ] = struct.type === 'function' ?
				methodScoping( scope, struct.val, struct ) :
				struct.type === 'object' ?
					{} :
					struct.val;
	
			obj[ struct.name ].__dt_wrapper = true;
	
			// Property extension
			_Api.extend( scope, obj[ struct.name ], struct.propExt );
		}
	};
	
	
	// @todo - Is there need for an augment function?
	// _Api.augment = function ( inst, name )
	// {
	// 	// Find src object in the structure from the name
	// 	var parts = name.split('.');
	
	// 	_Api.extend( inst, obj );
	// };
	
	
	//     [
	//       {
	//         name:      'data'                -- string   - Property name
	//         val:       function () {},       -- function - Api method (or undefined if just an object
	//         methodExt: [ ... ],              -- array    - Array of Api object definitions to extend the method result
	//         propExt:   [ ... ]               -- array    - Array of Api object definitions to extend the property
	//       },
	//       {
	//         name:     'row'
	//         val:       {},
	//         methodExt: [ ... ],
	//         propExt:   [
	//           {
	//             name:      'data'
	//             val:       function () {},
	//             methodExt: [ ... ],
	//             propExt:   [ ... ]
	//           },
	//           ...
	//         ]
	//       }
	//     ]
	
	_Api.register = _api_register = function ( name, val )
	{
		if ( Array.isArray( name ) ) {
			for ( var j=0, jen=name.length ; j<jen ; j++ ) {
				_Api.register( name[j], val );
			}
			return;
		}
	
		var
			i, ien,
			heir = name.split('.'),
			struct = __apiStruct,
			key, method;
	
		var find = function ( src, name ) {
			for ( var i=0, ien=src.length ; i<ien ; i++ ) {
				if ( src[i].name === name ) {
					return src[i];
				}
			}
			return null;
		};
	
		for ( i=0, ien=heir.length ; i<ien ; i++ ) {
			method = heir[i].indexOf('()') !== -1;
			key = method ?
				heir[i].replace('()', '') :
				heir[i];
	
			var src = find( struct, key );
			if ( ! src ) {
				src = {
					name:      key,
					val:       {},
					methodExt: [],
					propExt:   [],
					type:      'object'
				};
				struct.push( src );
			}
	
			if ( i === ien-1 ) {
				src.val = val;
				src.type = typeof val === 'function' ?
					'function' :
					$.isPlainObject( val ) ?
						'object' :
						'other';
			}
			else {
				struct = method ?
					src.methodExt :
					src.propExt;
			}
		}
	};
	
	_Api.registerPlural = _api_registerPlural = function ( pluralName, singularName, val ) {
		_Api.register( pluralName, val );
	
		_Api.register( singularName, function () {
			var ret = val.apply( this, arguments );
	
			if ( ret === this ) {
				// Returned item is the API instance that was passed in, return it
				return this;
			}
			else if ( ret instanceof _Api ) {
				// New API instance returned, want the value from the first item
				// in the returned array for the singular result.
				return ret.length ?
					Array.isArray( ret[0] ) ?
						new _Api( ret.context, ret[0] ) : // Array results are 'enhanced'
						ret[0] :
					undefined;
			}
	
			// Non-API return - just fire it back
			return ret;
		} );
	};
	
	
	/**
	 * Selector for HTML tables. Apply the given selector to the give array of
	 * DataTables settings objects.
	 *
	 * @param {string|integer} [selector] jQuery selector string or integer
	 * @param  {array} Array of DataTables settings objects to be filtered
	 * @return {array}
	 * @ignore
	 */
	var __table_selector = function ( selector, a )
	{
		if ( Array.isArray(selector) ) {
			return $.map( selector, function (item) {
				return __table_selector(item, a);
			} );
		}
	
		// Integer is used to pick out a table by index
		if ( typeof selector === 'number' ) {
			return [ a[ selector ] ];
		}
	
		// Perform a jQuery selector on the table nodes
		var nodes = $.map( a, function (el, i) {
			return el.nTable;
		} );
	
		return $(nodes)
			.filter( selector )
			.map( function (i) {
				// Need to translate back from the table node to the settings
				var idx = $.inArray( this, nodes );
				return a[ idx ];
			} )
			.toArray();
	};
	
	
	
	/**
	 * Context selector for the API's context (i.e. the tables the API instance
	 * refers to.
	 *
	 * @name    DataTable.Api#tables
	 * @param {string|integer} [selector] Selector to pick which tables the iterator
	 *   should operate on. If not given, all tables in the current context are
	 *   used. This can be given as a jQuery selector (for example `':gt(0)'`) to
	 *   select multiple tables or as an integer to select a single table.
	 * @returns {DataTable.Api} Returns a new API instance if a selector is given.
	 */
	_api_register( 'tables()', function ( selector ) {
		// A new instance is created if there was a selector specified
		return selector !== undefined && selector !== null ?
			new _Api( __table_selector( selector, this.context ) ) :
			this;
	} );
	
	
	_api_register( 'table()', function ( selector ) {
		var tables = this.tables( selector );
		var ctx = tables.context;
	
		// Truncate to the first matched table
		return ctx.length ?
			new _Api( ctx[0] ) :
			tables;
	} );
	
	
	_api_registerPlural( 'tables().nodes()', 'table().node()' , function () {
		return this.iterator( 'table', function ( ctx ) {
			return ctx.nTable;
		}, 1 );
	} );
	
	
	_api_registerPlural( 'tables().body()', 'table().body()' , function () {
		return this.iterator( 'table', function ( ctx ) {
			return ctx.nTBody;
		}, 1 );
	} );
	
	
	_api_registerPlural( 'tables().header()', 'table().header()' , function () {
		return this.iterator( 'table', function ( ctx ) {
			return ctx.nTHead;
		}, 1 );
	} );
	
	
	_api_registerPlural( 'tables().footer()', 'table().footer()' , function () {
		return this.iterator( 'table', function ( ctx ) {
			return ctx.nTFoot;
		}, 1 );
	} );
	
	
	_api_registerPlural( 'tables().containers()', 'table().container()' , function () {
		return this.iterator( 'table', function ( ctx ) {
			return ctx.nTableWrapper;
		}, 1 );
	} );
	
	
	
	/**
	 * Redraw the tables in the current context.
	 */
	_api_register( 'draw()', function ( paging ) {
		return this.iterator( 'table', function ( settings ) {
			if ( paging === 'page' ) {
				_fnDraw( settings );
			}
			else {
				if ( typeof paging === 'string' ) {
					paging = paging === 'full-hold' ?
						false :
						true;
				}
	
				_fnReDraw( settings, paging===false );
			}
		} );
	} );
	
	
	
	/**
	 * Get the current page index.
	 *
	 * @return {integer} Current page index (zero based)
	 *//**
	 * Set the current page.
	 *
	 * Note that if you attempt to show a page which does not exist, DataTables will
	 * not throw an error, but rather reset the paging.
	 *
	 * @param {integer|string} action The paging action to take. This can be one of:
	 *  * `integer` - The page index to jump to
	 *  * `string` - An action to take:
	 *    * `first` - Jump to first page.
	 *    * `next` - Jump to the next page
	 *    * `previous` - Jump to previous page
	 *    * `last` - Jump to the last page.
	 * @returns {DataTables.Api} this
	 */
	_api_register( 'page()', function ( action ) {
		if ( action === undefined ) {
			return this.page.info().page; // not an expensive call
		}
	
		// else, have an action to take on all tables
		return this.iterator( 'table', function ( settings ) {
			_fnPageChange( settings, action );
		} );
	} );
	
	
	/**
	 * Paging information for the first table in the current context.
	 *
	 * If you require paging information for another table, use the `table()` method
	 * with a suitable selector.
	 *
	 * @return {object} Object with the following properties set:
	 *  * `page` - Current page index (zero based - i.e. the first page is `0`)
	 *  * `pages` - Total number of pages
	 *  * `start` - Display index for the first record shown on the current page
	 *  * `end` - Display index for the last record shown on the current page
	 *  * `length` - Display length (number of records). Note that generally `start
	 *    + length = end`, but this is not always true, for example if there are
	 *    only 2 records to show on the final page, with a length of 10.
	 *  * `recordsTotal` - Full data set length
	 *  * `recordsDisplay` - Data set length once the current filtering criterion
	 *    are applied.
	 */
	_api_register( 'page.info()', function ( action ) {
		if ( this.context.length === 0 ) {
			return undefined;
		}
	
		var
			settings   = this.context[0],
			start      = settings._iDisplayStart,
			len        = settings.oFeatures.bPaginate ? settings._iDisplayLength : -1,
			visRecords = settings.fnRecordsDisplay(),
			all        = len === -1;
	
		return {
			"page":           all ? 0 : Math.floor( start / len ),
			"pages":          all ? 1 : Math.ceil( visRecords / len ),
			"start":          start,
			"end":            settings.fnDisplayEnd(),
			"length":         len,
			"recordsTotal":   settings.fnRecordsTotal(),
			"recordsDisplay": visRecords,
			"serverSide":     _fnDataSource( settings ) === 'ssp'
		};
	} );
	
	
	/**
	 * Get the current page length.
	 *
	 * @return {integer} Current page length. Note `-1` indicates that all records
	 *   are to be shown.
	 *//**
	 * Set the current page length.
	 *
	 * @param {integer} Page length to set. Use `-1` to show all records.
	 * @returns {DataTables.Api} this
	 */
	_api_register( 'page.len()', function ( len ) {
		// Note that we can't call this function 'length()' because `length`
		// is a Javascript property of functions which defines how many arguments
		// the function expects.
		if ( len === undefined ) {
			return this.context.length !== 0 ?
				this.context[0]._iDisplayLength :
				undefined;
		}
	
		// else, set the page length
		return this.iterator( 'table', function ( settings ) {
			_fnLengthChange( settings, len );
		} );
	} );
	
	
	
	var __reload = function ( settings, holdPosition, callback ) {
		// Use the draw event to trigger a callback
		if ( callback ) {
			var api = new _Api( settings );
	
			api.one( 'draw', function () {
				callback( api.ajax.json() );
			} );
		}
	
		if ( _fnDataSource( settings ) == 'ssp' ) {
			_fnReDraw( settings, holdPosition );
		}
		else {
			_fnProcessingDisplay( settings, true );
	
			// Cancel an existing request
			var xhr = settings.jqXHR;
			if ( xhr && xhr.readyState !== 4 ) {
				xhr.abort();
			}
	
			// Trigger xhr
			_fnBuildAjax( settings, [], function( json ) {
				_fnClearTable( settings );
	
				var data = _fnAjaxDataSrc( settings, json );
				for ( var i=0, ien=data.length ; i<ien ; i++ ) {
					_fnAddData( settings, data[i] );
				}
	
				_fnReDraw( settings, holdPosition );
				_fnProcessingDisplay( settings, false );
			} );
		}
	};
	
	
	/**
	 * Get the JSON response from the last Ajax request that DataTables made to the
	 * server. Note that this returns the JSON from the first table in the current
	 * context.
	 *
	 * @return {object} JSON received from the server.
	 */
	_api_register( 'ajax.json()', function () {
		var ctx = this.context;
	
		if ( ctx.length > 0 ) {
			return ctx[0].json;
		}
	
		// else return undefined;
	} );
	
	
	/**
	 * Get the data submitted in the last Ajax request
	 */
	_api_register( 'ajax.params()', function () {
		var ctx = this.context;
	
		if ( ctx.length > 0 ) {
			return ctx[0].oAjaxData;
		}
	
		// else return undefined;
	} );
	
	
	/**
	 * Reload tables from the Ajax data source. Note that this function will
	 * automatically re-draw the table when the remote data has been loaded.
	 *
	 * @param {boolean} [reset=true] Reset (default) or hold the current paging
	 *   position. A full re-sort and re-filter is performed when this method is
	 *   called, which is why the pagination reset is the default action.
	 * @returns {DataTables.Api} this
	 */
	_api_register( 'ajax.reload()', function ( callback, resetPaging ) {
		return this.iterator( 'table', function (settings) {
			__reload( settings, resetPaging===false, callback );
		} );
	} );
	
	
	/**
	 * Get the current Ajax URL. Note that this returns the URL from the first
	 * table in the current context.
	 *
	 * @return {string} Current Ajax source URL
	 *//**
	 * Set the Ajax URL. Note that this will set the URL for all tables in the
	 * current context.
	 *
	 * @param {string} url URL to set.
	 * @returns {DataTables.Api} this
	 */
	_api_register( 'ajax.url()', function ( url ) {
		var ctx = this.context;
	
		if ( url === undefined ) {
			// get
			if ( ctx.length === 0 ) {
				return undefined;
			}
			ctx = ctx[0];
	
			return ctx.ajax ?
				$.isPlainObject( ctx.ajax ) ?
					ctx.ajax.url :
					ctx.ajax :
				ctx.sAjaxSource;
		}
	
		// set
		return this.iterator( 'table', function ( settings ) {
			if ( $.isPlainObject( settings.ajax ) ) {
				settings.ajax.url = url;
			}
			else {
				settings.ajax = url;
			}
			// No need to consider sAjaxSource here since DataTables gives priority
			// to `ajax` over `sAjaxSource`. So setting `ajax` here, renders any
			// value of `sAjaxSource` redundant.
		} );
	} );
	
	
	/**
	 * Load data from the newly set Ajax URL. Note that this method is only
	 * available when `ajax.url()` is used to set a URL. Additionally, this method
	 * has the same effect as calling `ajax.reload()` but is provided for
	 * convenience when setting a new URL. Like `ajax.reload()` it will
	 * automatically redraw the table once the remote data has been loaded.
	 *
	 * @returns {DataTables.Api} this
	 */
	_api_register( 'ajax.url().load()', function ( callback, resetPaging ) {
		// Same as a reload, but makes sense to present it for easy access after a
		// url change
		return this.iterator( 'table', function ( ctx ) {
			__reload( ctx, resetPaging===false, callback );
		} );
	} );
	
	
	
	
	var _selector_run = function ( type, selector, selectFn, settings, opts )
	{
		var
			out = [], res,
			a, i, ien, j, jen,
			selectorType = typeof selector;
	
		// Can't just check for isArray here, as an API or jQuery instance might be
		// given with their array like look
		if ( ! selector || selectorType === 'string' || selectorType === 'function' || selector.length === undefined ) {
			selector = [ selector ];
		}
	
		for ( i=0, ien=selector.length ; i<ien ; i++ ) {
			// Only split on simple strings - complex expressions will be jQuery selectors
			a = selector[i] && selector[i].split && ! selector[i].match(/[\[\(:]/) ?
				selector[i].split(',') :
				[ selector[i] ];
	
			for ( j=0, jen=a.length ; j<jen ; j++ ) {
				res = selectFn( typeof a[j] === 'string' ? (a[j]).trim() : a[j] );
	
				if ( res && res.length ) {
					out = out.concat( res );
				}
			}
		}
	
		// selector extensions
		var ext = _ext.selector[ type ];
		if ( ext.length ) {
			for ( i=0, ien=ext.length ; i<ien ; i++ ) {
				out = ext[i]( settings, opts, out );
			}
		}
	
		return _unique( out );
	};
	
	
	var _selector_opts = function ( opts )
	{
		if ( ! opts ) {
			opts = {};
		}
	
		// Backwards compatibility for 1.9- which used the terminology filter rather
		// than search
		if ( opts.filter && opts.search === undefined ) {
			opts.search = opts.filter;
		}
	
		return $.extend( {
			search: 'none',
			order: 'current',
			page: 'all'
		}, opts );
	};
	
	
	var _selector_first = function ( inst )
	{
		// Reduce the API instance to the first item found
		for ( var i=0, ien=inst.length ; i<ien ; i++ ) {
			if ( inst[i].length > 0 ) {
				// Assign the first element to the first item in the instance
				// and truncate the instance and context
				inst[0] = inst[i];
				inst[0].length = 1;
				inst.length = 1;
				inst.context = [ inst.context[i] ];
	
				return inst;
			}
		}
	
		// Not found - return an empty instance
		inst.length = 0;
		return inst;
	};
	
	
	var _selector_row_indexes = function ( settings, opts )
	{
		var
			i, ien, tmp, a=[],
			displayFiltered = settings.aiDisplay,
			displayMaster = settings.aiDisplayMaster;
	
		var
			search = opts.search,  // none, applied, removed
			order  = opts.order,   // applied, current, index (original - compatibility with 1.9)
			page   = opts.page;    // all, current
	
		if ( _fnDataSource( settings ) == 'ssp' ) {
			// In server-side processing mode, most options are irrelevant since
			// rows not shown don't exist and the index order is the applied order
			// Removed is a special case - for consistency just return an empty
			// array
			return search === 'removed' ?
				[] :
				_range( 0, displayMaster.length );
		}
		else if ( page == 'current' ) {
			// Current page implies that order=current and filter=applied, since it is
			// fairly senseless otherwise, regardless of what order and search actually
			// are
			for ( i=settings._iDisplayStart, ien=settings.fnDisplayEnd() ; i<ien ; i++ ) {
				a.push( displayFiltered[i] );
			}
		}
		else if ( order == 'current' || order == 'applied' ) {
			if ( search == 'none') {
				a = displayMaster.slice();
			}
			else if ( search == 'applied' ) {
				a = displayFiltered.slice();
			}
			else if ( search == 'removed' ) {
				// O(n+m) solution by creating a hash map
				var displayFilteredMap = {};
	
				for ( var i=0, ien=displayFiltered.length ; i<ien ; i++ ) {
					displayFilteredMap[displayFiltered[i]] = null;
				}
	
				a = $.map( displayMaster, function (el) {
					return ! displayFilteredMap.hasOwnProperty(el) ?
						el :
						null;
				} );
			}
		}
		else if ( order == 'index' || order == 'original' ) {
			for ( i=0, ien=settings.aoData.length ; i<ien ; i++ ) {
				if ( search == 'none' ) {
					a.push( i );
				}
				else { // applied | removed
					tmp = $.inArray( i, displayFiltered );
	
					if ((tmp === -1 && search == 'removed') ||
						(tmp >= 0   && search == 'applied') )
					{
						a.push( i );
					}
				}
			}
		}
	
		return a;
	};
	
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Rows
	 *
	 * {}          - no selector - use all available rows
	 * {integer}   - row aoData index
	 * {node}      - TR node
	 * {string}    - jQuery selector to apply to the TR elements
	 * {array}     - jQuery array of nodes, or simply an array of TR nodes
	 *
	 */
	var __row_selector = function ( settings, selector, opts )
	{
		var rows;
		var run = function ( sel ) {
			var selInt = _intVal( sel );
			var i, ien;
			var aoData = settings.aoData;
	
			// Short cut - selector is a number and no options provided (default is
			// all records, so no need to check if the index is in there, since it
			// must be - dev error if the index doesn't exist).
			if ( selInt !== null && ! opts ) {
				return [ selInt ];
			}
	
			if ( ! rows ) {
				rows = _selector_row_indexes( settings, opts );
			}
	
			if ( selInt !== null && $.inArray( selInt, rows ) !== -1 ) {
				// Selector - integer
				return [ selInt ];
			}
			else if ( sel === null || sel === undefined || sel === '' ) {
				// Selector - none
				return rows;
			}
	
			// Selector - function
			if ( typeof sel === 'function' ) {
				return $.map( rows, function (idx) {
					var row = aoData[ idx ];
					return sel( idx, row._aData, row.nTr ) ? idx : null;
				} );
			}
	
			// Selector - node
			if ( sel.nodeName ) {
				var rowIdx = sel._DT_RowIndex;  // Property added by DT for fast lookup
				var cellIdx = sel._DT_CellIndex;
	
				if ( rowIdx !== undefined ) {
					// Make sure that the row is actually still present in the table
					return aoData[ rowIdx ] && aoData[ rowIdx ].nTr === sel ?
						[ rowIdx ] :
						[];
				}
				else if ( cellIdx ) {
					return aoData[ cellIdx.row ] && aoData[ cellIdx.row ].nTr === sel.parentNode ?
						[ cellIdx.row ] :
						[];
				}
				else {
					var host = $(sel).closest('*[data-dt-row]');
					return host.length ?
						[ host.data('dt-row') ] :
						[];
				}
			}
	
			// ID selector. Want to always be able to select rows by id, regardless
			// of if the tr element has been created or not, so can't rely upon
			// jQuery here - hence a custom implementation. This does not match
			// Sizzle's fast selector or HTML4 - in HTML5 the ID can be anything,
			// but to select it using a CSS selector engine (like Sizzle or
			// querySelect) it would need to need to be escaped for some characters.
			// DataTables simplifies this for row selectors since you can select
			// only a row. A # indicates an id any anything that follows is the id -
			// unescaped.
			if ( typeof sel === 'string' && sel.charAt(0) === '#' ) {
				// get row index from id
				var rowObj = settings.aIds[ sel.replace( /^#/, '' ) ];
				if ( rowObj !== undefined ) {
					return [ rowObj.idx ];
				}
	
				// need to fall through to jQuery in case there is DOM id that
				// matches
			}
			
			// Get nodes in the order from the `rows` array with null values removed
			var nodes = _removeEmpty(
				_pluck_order( settings.aoData, rows, 'nTr' )
			);
	
			// Selector - jQuery selector string, array of nodes or jQuery object/
			// As jQuery's .filter() allows jQuery objects to be passed in filter,
			// it also allows arrays, so this will cope with all three options
			return $(nodes)
				.filter( sel )
				.map( function () {
					return this._DT_RowIndex;
				} )
				.toArray();
		};
	
		return _selector_run( 'row', selector, run, settings, opts );
	};
	
	
	_api_register( 'rows()', function ( selector, opts ) {
		// argument shifting
		if ( selector === undefined ) {
			selector = '';
		}
		else if ( $.isPlainObject( selector ) ) {
			opts = selector;
			selector = '';
		}
	
		opts = _selector_opts( opts );
	
		var inst = this.iterator( 'table', function ( settings ) {
			return __row_selector( settings, selector, opts );
		}, 1 );
	
		// Want argument shifting here and in __row_selector?
		inst.selector.rows = selector;
		inst.selector.opts = opts;
	
		return inst;
	} );
	
	_api_register( 'rows().nodes()', function () {
		return this.iterator( 'row', function ( settings, row ) {
			return settings.aoData[ row ].nTr || undefined;
		}, 1 );
	} );
	
	_api_register( 'rows().data()', function () {
		return this.iterator( true, 'rows', function ( settings, rows ) {
			return _pluck_order( settings.aoData, rows, '_aData' );
		}, 1 );
	} );
	
	_api_registerPlural( 'rows().cache()', 'row().cache()', function ( type ) {
		return this.iterator( 'row', function ( settings, row ) {
			var r = settings.aoData[ row ];
			return type === 'search' ? r._aFilterData : r._aSortData;
		}, 1 );
	} );
	
	_api_registerPlural( 'rows().invalidate()', 'row().invalidate()', function ( src ) {
		return this.iterator( 'row', function ( settings, row ) {
			_fnInvalidate( settings, row, src );
		} );
	} );
	
	_api_registerPlural( 'rows().indexes()', 'row().index()', function () {
		return this.iterator( 'row', function ( settings, row ) {
			return row;
		}, 1 );
	} );
	
	_api_registerPlural( 'rows().ids()', 'row().id()', function ( hash ) {
		var a = [];
		var context = this.context;
	
		// `iterator` will drop undefined values, but in this case we want them
		for ( var i=0, ien=context.length ; i<ien ; i++ ) {
			for ( var j=0, jen=this[i].length ; j<jen ; j++ ) {
				var id = context[i].rowIdFn( context[i].aoData[ this[i][j] ]._aData );
				a.push( (hash === true ? '#' : '' )+ id );
			}
		}
	
		return new _Api( context, a );
	} );
	
	_api_registerPlural( 'rows().remove()', 'row().remove()', function () {
		var that = this;
	
		this.iterator( 'row', function ( settings, row, thatIdx ) {
			var data = settings.aoData;
			var rowData = data[ row ];
			var i, ien, j, jen;
			var loopRow, loopCells;
	
			data.splice( row, 1 );
	
			// Update the cached indexes
			for ( i=0, ien=data.length ; i<ien ; i++ ) {
				loopRow = data[i];
				loopCells = loopRow.anCells;
	
				// Rows
				if ( loopRow.nTr !== null ) {
					loopRow.nTr._DT_RowIndex = i;
				}
	
				// Cells
				if ( loopCells !== null ) {
					for ( j=0, jen=loopCells.length ; j<jen ; j++ ) {
						loopCells[j]._DT_CellIndex.row = i;
					}
				}
			}
	
			// Delete from the display arrays
			_fnDeleteIndex( settings.aiDisplayMaster, row );
			_fnDeleteIndex( settings.aiDisplay, row );
			_fnDeleteIndex( that[ thatIdx ], row, false ); // maintain local indexes
	
			// For server-side processing tables - subtract the deleted row from the count
			if ( settings._iRecordsDisplay > 0 ) {
				settings._iRecordsDisplay--;
			}
	
			// Check for an 'overflow' they case for displaying the table
			_fnLengthOverflow( settings );
	
			// Remove the row's ID reference if there is one
			var id = settings.rowIdFn( rowData._aData );
			if ( id !== undefined ) {
				delete settings.aIds[ id ];
			}
		} );
	
		this.iterator( 'table', function ( settings ) {
			for ( var i=0, ien=settings.aoData.length ; i<ien ; i++ ) {
				settings.aoData[i].idx = i;
			}
		} );
	
		return this;
	} );
	
	
	_api_register( 'rows.add()', function ( rows ) {
		var newRows = this.iterator( 'table', function ( settings ) {
				var row, i, ien;
				var out = [];
	
				for ( i=0, ien=rows.length ; i<ien ; i++ ) {
					row = rows[i];
	
					if ( row.nodeName && row.nodeName.toUpperCase() === 'TR' ) {
						out.push( _fnAddTr( settings, row )[0] );
					}
					else {
						out.push( _fnAddData( settings, row ) );
					}
				}
	
				return out;
			}, 1 );
	
		// Return an Api.rows() extended instance, so rows().nodes() etc can be used
		var modRows = this.rows( -1 );
		modRows.pop();
		$.merge( modRows, newRows );
	
		return modRows;
	} );
	
	
	
	
	
	/**
	 *
	 */
	_api_register( 'row()', function ( selector, opts ) {
		return _selector_first( this.rows( selector, opts ) );
	} );
	
	
	_api_register( 'row().data()', function ( data ) {
		var ctx = this.context;
	
		if ( data === undefined ) {
			// Get
			return ctx.length && this.length ?
				ctx[0].aoData[ this[0] ]._aData :
				undefined;
		}
	
		// Set
		var row = ctx[0].aoData[ this[0] ];
		row._aData = data;
	
		// If the DOM has an id, and the data source is an array
		if ( Array.isArray( data ) && row.nTr && row.nTr.id ) {
			_fnSetObjectDataFn( ctx[0].rowId )( data, row.nTr.id );
		}
	
		// Automatically invalidate
		_fnInvalidate( ctx[0], this[0], 'data' );
	
		return this;
	} );
	
	
	_api_register( 'row().node()', function () {
		var ctx = this.context;
	
		return ctx.length && this.length ?
			ctx[0].aoData[ this[0] ].nTr || null :
			null;
	} );
	
	
	_api_register( 'row.add()', function ( row ) {
		// Allow a jQuery object to be passed in - only a single row is added from
		// it though - the first element in the set
		if ( row instanceof $ && row.length ) {
			row = row[0];
		}
	
		var rows = this.iterator( 'table', function ( settings ) {
			if ( row.nodeName && row.nodeName.toUpperCase() === 'TR' ) {
				return _fnAddTr( settings, row )[0];
			}
			return _fnAddData( settings, row );
		} );
	
		// Return an Api.rows() extended instance, with the newly added row selected
		return this.row( rows[0] );
	} );
	
	
	$(document).on('plugin-init.dt', function (e, context) {
		var api = new _Api( context );
		var namespace = 'on-plugin-init';
		var stateSaveParamsEvent = 'stateSaveParams.' + namespace;
		var destroyEvent = 'destroy. ' + namespace;
	
		api.on( stateSaveParamsEvent, function ( e, settings, d ) {
			// This could be more compact with the API, but it is a lot faster as a simple
			// internal loop
			var idFn = settings.rowIdFn;
			var data = settings.aoData;
			var ids = [];
	
			for (var i=0 ; i<data.length ; i++) {
				if (data[i]._detailsShow) {
					ids.push( '#' + idFn(data[i]._aData) );
				}
			}
	
			d.childRows = ids;
		});
	
		api.on( destroyEvent, function () {
			api.off(stateSaveParamsEvent + ' ' + destroyEvent);
		});
	
		var loaded = api.state.loaded();
	
		if ( loaded && loaded.childRows ) {
			api
				.rows( $.map(loaded.childRows, function (id){
					return id.replace(/:/g, '\\:')
				}) )
				.every( function () {
					_fnCallbackFire( context, null, 'requestChild', [ this ] )
				});
		}
	});
	
	var __details_add = function ( ctx, row, data, klass )
	{
		// Convert to array of TR elements
		var rows = [];
		var addRow = function ( r, k ) {
			// Recursion to allow for arrays of jQuery objects
			if ( Array.isArray( r ) || r instanceof $ ) {
				for ( var i=0, ien=r.length ; i<ien ; i++ ) {
					addRow( r[i], k );
				}
				return;
			}
	
			// If we get a TR element, then just add it directly - up to the dev
			// to add the correct number of columns etc
			if ( r.nodeName && r.nodeName.toLowerCase() === 'tr' ) {
				rows.push( r );
			}
			else {
				// Otherwise create a row with a wrapper
				var created = $('<tr><td></td></tr>').addClass( k );
				$('td', created)
					.addClass( k )
					.html( r )
					[0].colSpan = _fnVisbleColumns( ctx );
	
				rows.push( created[0] );
			}
		};
	
		addRow( data, klass );
	
		if ( row._details ) {
			row._details.detach();
		}
	
		row._details = $(rows);
	
		// If the children were already shown, that state should be retained
		if ( row._detailsShow ) {
			row._details.insertAfter( row.nTr );
		}
	};
	
	
	// Make state saving of child row details async to allow them to be batch processed
	var __details_state = DataTable.util.throttle(
		function (ctx) {
			_fnSaveState( ctx[0] )
		},
		500
	);
	
	
	var __details_remove = function ( api, idx )
	{
		var ctx = api.context;
	
		if ( ctx.length ) {
			var row = ctx[0].aoData[ idx !== undefined ? idx : api[0] ];
	
			if ( row && row._details ) {
				row._details.remove();
	
				row._detailsShow = undefined;
				row._details = undefined;
				$( row.nTr ).removeClass( 'dt-hasChild' );
				__details_state( ctx );
			}
		}
	};
	
	
	var __details_display = function ( api, show ) {
		var ctx = api.context;
	
		if ( ctx.length && api.length ) {
			var row = ctx[0].aoData[ api[0] ];
	
			if ( row._details ) {
				row._detailsShow = show;
	
				if ( show ) {
					row._details.insertAfter( row.nTr );
					$( row.nTr ).addClass( 'dt-hasChild' );
				}
				else {
					row._details.detach();
					$( row.nTr ).removeClass( 'dt-hasChild' );
				}
	
				_fnCallbackFire( ctx[0], null, 'childRow', [ show, api.row( api[0] ) ] )
	
				__details_events( ctx[0] );
				__details_state( ctx );
			}
		}
	};
	
	
	var __details_events = function ( settings )
	{
		var api = new _Api( settings );
		var namespace = '.dt.DT_details';
		var drawEvent = 'draw'+namespace;
		var colvisEvent = 'column-sizing'+namespace;
		var destroyEvent = 'destroy'+namespace;
		var data = settings.aoData;
	
		api.off( drawEvent +' '+ colvisEvent +' '+ destroyEvent );
	
		if ( _pluck( data, '_details' ).length > 0 ) {
			// On each draw, insert the required elements into the document
			api.on( drawEvent, function ( e, ctx ) {
				if ( settings !== ctx ) {
					return;
				}
	
				api.rows( {page:'current'} ).eq(0).each( function (idx) {
					// Internal data grab
					var row = data[ idx ];
	
					if ( row._detailsShow ) {
						row._details.insertAfter( row.nTr );
					}
				} );
			} );
	
			// Column visibility change - update the colspan
			api.on( colvisEvent, function ( e, ctx, idx, vis ) {
				if ( settings !== ctx ) {
					return;
				}
	
				// Update the colspan for the details rows (note, only if it already has
				// a colspan)
				var row, visible = _fnVisbleColumns( ctx );
	
				for ( var i=0, ien=data.length ; i<ien ; i++ ) {
					row = data[i];
	
					if ( row._details ) {
						row._details.children('td[colspan]').attr('colspan', visible );
					}
				}
			} );
	
			// Table destroyed - nuke any child rows
			api.on( destroyEvent, function ( e, ctx ) {
				if ( settings !== ctx ) {
					return;
				}
	
				for ( var i=0, ien=data.length ; i<ien ; i++ ) {
					if ( data[i]._details ) {
						__details_remove( api, i );
					}
				}
			} );
		}
	};
	
	// Strings for the method names to help minification
	var _emp = '';
	var _child_obj = _emp+'row().child';
	var _child_mth = _child_obj+'()';
	
	// data can be:
	//  tr
	//  string
	//  jQuery or array of any of the above
	_api_register( _child_mth, function ( data, klass ) {
		var ctx = this.context;
	
		if ( data === undefined ) {
			// get
			return ctx.length && this.length ?
				ctx[0].aoData[ this[0] ]._details :
				undefined;
		}
		else if ( data === true ) {
			// show
			this.child.show();
		}
		else if ( data === false ) {
			// remove
			__details_remove( this );
		}
		else if ( ctx.length && this.length ) {
			// set
			__details_add( ctx[0], ctx[0].aoData[ this[0] ], data, klass );
		}
	
		return this;
	} );
	
	
	_api_register( [
		_child_obj+'.show()',
		_child_mth+'.show()' // only when `child()` was called with parameters (without
	], function ( show ) {   // it returns an object and this method is not executed)
		__details_display( this, true );
		return this;
	} );
	
	
	_api_register( [
		_child_obj+'.hide()',
		_child_mth+'.hide()' // only when `child()` was called with parameters (without
	], function () {         // it returns an object and this method is not executed)
		__details_display( this, false );
		return this;
	} );
	
	
	_api_register( [
		_child_obj+'.remove()',
		_child_mth+'.remove()' // only when `child()` was called with parameters (without
	], function () {           // it returns an object and this method is not executed)
		__details_remove( this );
		return this;
	} );
	
	
	_api_register( _child_obj+'.isShown()', function () {
		var ctx = this.context;
	
		if ( ctx.length && this.length ) {
			// _detailsShown as false or undefined will fall through to return false
			return ctx[0].aoData[ this[0] ]._detailsShow || false;
		}
		return false;
	} );
	
	
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Columns
	 *
	 * {integer}           - column index (>=0 count from left, <0 count from right)
	 * "{integer}:visIdx"  - visible column index (i.e. translate to column index)  (>=0 count from left, <0 count from right)
	 * "{integer}:visible" - alias for {integer}:visIdx  (>=0 count from left, <0 count from right)
	 * "{string}:name"     - column name
	 * "{string}"          - jQuery selector on column header nodes
	 *
	 */
	
	// can be an array of these items, comma separated list, or an array of comma
	// separated lists
	
	var __re_column_selector = /^([^:]+):(name|visIdx|visible)$/;
	
	
	// r1 and r2 are redundant - but it means that the parameters match for the
	// iterator callback in columns().data()
	var __columnData = function ( settings, column, r1, r2, rows ) {
		var a = [];
		for ( var row=0, ien=rows.length ; row<ien ; row++ ) {
			a.push( _fnGetCellData( settings, rows[row], column ) );
		}
		return a;
	};
	
	
	var __column_selector = function ( settings, selector, opts )
	{
		var
			columns = settings.aoColumns,
			names = _pluck( columns, 'sName' ),
			nodes = _pluck( columns, 'nTh' );
	
		var run = function ( s ) {
			var selInt = _intVal( s );
	
			// Selector - all
			if ( s === '' ) {
				return _range( columns.length );
			}
	
			// Selector - index
			if ( selInt !== null ) {
				return [ selInt >= 0 ?
					selInt : // Count from left
					columns.length + selInt // Count from right (+ because its a negative value)
				];
			}
	
			// Selector = function
			if ( typeof s === 'function' ) {
				var rows = _selector_row_indexes( settings, opts );
	
				return $.map( columns, function (col, idx) {
					return s(
							idx,
							__columnData( settings, idx, 0, 0, rows ),
							nodes[ idx ]
						) ? idx : null;
				} );
			}
	
			// jQuery or string selector
			var match = typeof s === 'string' ?
				s.match( __re_column_selector ) :
				'';
	
			if ( match ) {
				switch( match[2] ) {
					case 'visIdx':
					case 'visible':
						var idx = parseInt( match[1], 10 );
						// Visible index given, convert to column index
						if ( idx < 0 ) {
							// Counting from the right
							var visColumns = $.map( columns, function (col,i) {
								return col.bVisible ? i : null;
							} );
							return [ visColumns[ visColumns.length + idx ] ];
						}
						// Counting from the left
						return [ _fnVisibleToColumnIndex( settings, idx ) ];
	
					case 'name':
						// match by name. `names` is column index complete and in order
						return $.map( names, function (name, i) {
							return name === match[1] ? i : null;
						} );
	
					default:
						return [];
				}
			}
	
			// Cell in the table body
			if ( s.nodeName && s._DT_CellIndex ) {
				return [ s._DT_CellIndex.column ];
			}
	
			// jQuery selector on the TH elements for the columns
			var jqResult = $( nodes )
				.filter( s )
				.map( function () {
					return $.inArray( this, nodes ); // `nodes` is column index complete and in order
				} )
				.toArray();
	
			if ( jqResult.length || ! s.nodeName ) {
				return jqResult;
			}
	
			// Otherwise a node which might have a `dt-column` data attribute, or be
			// a child or such an element
			var host = $(s).closest('*[data-dt-column]');
			return host.length ?
				[ host.data('dt-column') ] :
				[];
		};
	
		return _selector_run( 'column', selector, run, settings, opts );
	};
	
	
	var __setColumnVis = function ( settings, column, vis ) {
		var
			cols = settings.aoColumns,
			col  = cols[ column ],
			data = settings.aoData,
			row, cells, i, ien, tr;
	
		// Get
		if ( vis === undefined ) {
			return col.bVisible;
		}
	
		// Set
		// No change
		if ( col.bVisible === vis ) {
			return;
		}
	
		if ( vis ) {
			// Insert column
			// Need to decide if we should use appendChild or insertBefore
			var insertBefore = $.inArray( true, _pluck(cols, 'bVisible'), column+1 );
	
			for ( i=0, ien=data.length ; i<ien ; i++ ) {
				tr = data[i].nTr;
				cells = data[i].anCells;
	
				if ( tr ) {
					// insertBefore can act like appendChild if 2nd arg is null
					tr.insertBefore( cells[ column ], cells[ insertBefore ] || null );
				}
			}
		}
		else {
			// Remove column
			$( _pluck( settings.aoData, 'anCells', column ) ).detach();
		}
	
		// Common actions
		col.bVisible = vis;
	};
	
	
	_api_register( 'columns()', function ( selector, opts ) {
		// argument shifting
		if ( selector === undefined ) {
			selector = '';
		}
		else if ( $.isPlainObject( selector ) ) {
			opts = selector;
			selector = '';
		}
	
		opts = _selector_opts( opts );
	
		var inst = this.iterator( 'table', function ( settings ) {
			return __column_selector( settings, selector, opts );
		}, 1 );
	
		// Want argument shifting here and in _row_selector?
		inst.selector.cols = selector;
		inst.selector.opts = opts;
	
		return inst;
	} );
	
	_api_registerPlural( 'columns().header()', 'column().header()', function ( selector, opts ) {
		return this.iterator( 'column', function ( settings, column ) {
			return settings.aoColumns[column].nTh;
		}, 1 );
	} );
	
	_api_registerPlural( 'columns().footer()', 'column().footer()', function ( selector, opts ) {
		return this.iterator( 'column', function ( settings, column ) {
			return settings.aoColumns[column].nTf;
		}, 1 );
	} );
	
	_api_registerPlural( 'columns().data()', 'column().data()', function () {
		return this.iterator( 'column-rows', __columnData, 1 );
	} );
	
	_api_registerPlural( 'columns().dataSrc()', 'column().dataSrc()', function () {
		return this.iterator( 'column', function ( settings, column ) {
			return settings.aoColumns[column].mData;
		}, 1 );
	} );
	
	_api_registerPlural( 'columns().cache()', 'column().cache()', function ( type ) {
		return this.iterator( 'column-rows', function ( settings, column, i, j, rows ) {
			return _pluck_order( settings.aoData, rows,
				type === 'search' ? '_aFilterData' : '_aSortData', column
			);
		}, 1 );
	} );
	
	_api_registerPlural( 'columns().nodes()', 'column().nodes()', function () {
		return this.iterator( 'column-rows', function ( settings, column, i, j, rows ) {
			return _pluck_order( settings.aoData, rows, 'anCells', column ) ;
		}, 1 );
	} );
	
	_api_registerPlural( 'columns().visible()', 'column().visible()', function ( vis, calc ) {
		var that = this;
		var ret = this.iterator( 'column', function ( settings, column ) {
			if ( vis === undefined ) {
				return settings.aoColumns[ column ].bVisible;
			} // else
			__setColumnVis( settings, column, vis );
		} );
	
		// Group the column visibility changes
		if ( vis !== undefined ) {
			this.iterator( 'table', function ( settings ) {
				// Redraw the header after changes
				_fnDrawHead( settings, settings.aoHeader );
				_fnDrawHead( settings, settings.aoFooter );
		
				// Update colspan for no records display. Child rows and extensions will use their own
				// listeners to do this - only need to update the empty table item here
				if ( ! settings.aiDisplay.length ) {
					$(settings.nTBody).find('td[colspan]').attr('colspan', _fnVisbleColumns(settings));
				}
		
				_fnSaveState( settings );
	
				// Second loop once the first is done for events
				that.iterator( 'column', function ( settings, column ) {
					_fnCallbackFire( settings, null, 'column-visibility', [settings, column, vis, calc] );
				} );
	
				if ( calc === undefined || calc ) {
					that.columns.adjust();
				}
			});
		}
	
		return ret;
	} );
	
	_api_registerPlural( 'columns().indexes()', 'column().index()', function ( type ) {
		return this.iterator( 'column', function ( settings, column ) {
			return type === 'visible' ?
				_fnColumnIndexToVisible( settings, column ) :
				column;
		}, 1 );
	} );
	
	_api_register( 'columns.adjust()', function () {
		return this.iterator( 'table', function ( settings ) {
			_fnAdjustColumnSizing( settings );
		}, 1 );
	} );
	
	_api_register( 'column.index()', function ( type, idx ) {
		if ( this.context.length !== 0 ) {
			var ctx = this.context[0];
	
			if ( type === 'fromVisible' || type === 'toData' ) {
				return _fnVisibleToColumnIndex( ctx, idx );
			}
			else if ( type === 'fromData' || type === 'toVisible' ) {
				return _fnColumnIndexToVisible( ctx, idx );
			}
		}
	} );
	
	_api_register( 'column()', function ( selector, opts ) {
		return _selector_first( this.columns( selector, opts ) );
	} );
	
	var __cell_selector = function ( settings, selector, opts )
	{
		var data = settings.aoData;
		var rows = _selector_row_indexes( settings, opts );
		var cells = _removeEmpty( _pluck_order( data, rows, 'anCells' ) );
		var allCells = $(_flatten( [], cells ));
		var row;
		var columns = settings.aoColumns.length;
		var a, i, ien, j, o, host;
	
		var run = function ( s ) {
			var fnSelector = typeof s === 'function';
	
			if ( s === null || s === undefined || fnSelector ) {
				// All cells and function selectors
				a = [];
	
				for ( i=0, ien=rows.length ; i<ien ; i++ ) {
					row = rows[i];
	
					for ( j=0 ; j<columns ; j++ ) {
						o = {
							row: row,
							column: j
						};
	
						if ( fnSelector ) {
							// Selector - function
							host = data[ row ];
	
							if ( s( o, _fnGetCellData(settings, row, j), host.anCells ? host.anCells[j] : null ) ) {
								a.push( o );
							}
						}
						else {
							// Selector - all
							a.push( o );
						}
					}
				}
	
				return a;
			}
			
			// Selector - index
			if ( $.isPlainObject( s ) ) {
				// Valid cell index and its in the array of selectable rows
				return s.column !== undefined && s.row !== undefined && $.inArray( s.row, rows ) !== -1 ?
					[s] :
					[];
			}
	
			// Selector - jQuery filtered cells
			var jqResult = allCells
				.filter( s )
				.map( function (i, el) {
					return { // use a new object, in case someone changes the values
						row:    el._DT_CellIndex.row,
						column: el._DT_CellIndex.column
	 				};
				} )
				.toArray();
	
			if ( jqResult.length || ! s.nodeName ) {
				return jqResult;
			}
	
			// Otherwise the selector is a node, and there is one last option - the
			// element might be a child of an element which has dt-row and dt-column
			// data attributes
			host = $(s).closest('*[data-dt-row]');
			return host.length ?
				[ {
					row: host.data('dt-row'),
					column: host.data('dt-column')
				} ] :
				[];
		};
	
		return _selector_run( 'cell', selector, run, settings, opts );
	};
	
	
	
	
	_api_register( 'cells()', function ( rowSelector, columnSelector, opts ) {
		// Argument shifting
		if ( $.isPlainObject( rowSelector ) ) {
			// Indexes
			if ( rowSelector.row === undefined ) {
				// Selector options in first parameter
				opts = rowSelector;
				rowSelector = null;
			}
			else {
				// Cell index objects in first parameter
				opts = columnSelector;
				columnSelector = null;
			}
		}
		if ( $.isPlainObject( columnSelector ) ) {
			opts = columnSelector;
			columnSelector = null;
		}
	
		// Cell selector
		if ( columnSelector === null || columnSelector === undefined ) {
			return this.iterator( 'table', function ( settings ) {
				return __cell_selector( settings, rowSelector, _selector_opts( opts ) );
			} );
		}
	
		// The default built in options need to apply to row and columns
		var internalOpts = opts ? {
			page: opts.page,
			order: opts.order,
			search: opts.search
		} : {};
	
		// Row + column selector
		var columns = this.columns( columnSelector, internalOpts );
		var rows = this.rows( rowSelector, internalOpts );
		var i, ien, j, jen;
	
		var cellsNoOpts = this.iterator( 'table', function ( settings, idx ) {
			var a = [];
	
			for ( i=0, ien=rows[idx].length ; i<ien ; i++ ) {
				for ( j=0, jen=columns[idx].length ; j<jen ; j++ ) {
					a.push( {
						row:    rows[idx][i],
						column: columns[idx][j]
					} );
				}
			}
	
			return a;
		}, 1 );
	
		// There is currently only one extension which uses a cell selector extension
		// It is a _major_ performance drag to run this if it isn't needed, so this is
		// an extension specific check at the moment
		var cells = opts && opts.selected ?
			this.cells( cellsNoOpts, opts ) :
			cellsNoOpts;
	
		$.extend( cells.selector, {
			cols: columnSelector,
			rows: rowSelector,
			opts: opts
		} );
	
		return cells;
	} );
	
	
	_api_registerPlural( 'cells().nodes()', 'cell().node()', function () {
		return this.iterator( 'cell', function ( settings, row, column ) {
			var data = settings.aoData[ row ];
	
			return data && data.anCells ?
				data.anCells[ column ] :
				undefined;
		}, 1 );
	} );
	
	
	_api_register( 'cells().data()', function () {
		return this.iterator( 'cell', function ( settings, row, column ) {
			return _fnGetCellData( settings, row, column );
		}, 1 );
	} );
	
	
	_api_registerPlural( 'cells().cache()', 'cell().cache()', function ( type ) {
		type = type === 'search' ? '_aFilterData' : '_aSortData';
	
		return this.iterator( 'cell', function ( settings, row, column ) {
			return settings.aoData[ row ][ type ][ column ];
		}, 1 );
	} );
	
	
	_api_registerPlural( 'cells().render()', 'cell().render()', function ( type ) {
		return this.iterator( 'cell', function ( settings, row, column ) {
			return _fnGetCellData( settings, row, column, type );
		}, 1 );
	} );
	
	
	_api_registerPlural( 'cells().indexes()', 'cell().index()', function () {
		return this.iterator( 'cell', function ( settings, row, column ) {
			return {
				row: row,
				column: column,
				columnVisible: _fnColumnIndexToVisible( settings, column )
			};
		}, 1 );
	} );
	
	
	_api_registerPlural( 'cells().invalidate()', 'cell().invalidate()', function ( src ) {
		return this.iterator( 'cell', function ( settings, row, column ) {
			_fnInvalidate( settings, row, src, column );
		} );
	} );
	
	
	
	_api_register( 'cell()', function ( rowSelector, columnSelector, opts ) {
		return _selector_first( this.cells( rowSelector, columnSelector, opts ) );
	} );
	
	
	_api_register( 'cell().data()', function ( data ) {
		var ctx = this.context;
		var cell = this[0];
	
		if ( data === undefined ) {
			// Get
			return ctx.length && cell.length ?
				_fnGetCellData( ctx[0], cell[0].row, cell[0].column ) :
				undefined;
		}
	
		// Set
		_fnSetCellData( ctx[0], cell[0].row, cell[0].column, data );
		_fnInvalidate( ctx[0], cell[0].row, 'data', cell[0].column );
	
		return this;
	} );
	
	
	
	/**
	 * Get current ordering (sorting) that has been applied to the table.
	 *
	 * @returns {array} 2D array containing the sorting information for the first
	 *   table in the current context. Each element in the parent array represents
	 *   a column being sorted upon (i.e. multi-sorting with two columns would have
	 *   2 inner arrays). The inner arrays may have 2 or 3 elements. The first is
	 *   the column index that the sorting condition applies to, the second is the
	 *   direction of the sort (`desc` or `asc`) and, optionally, the third is the
	 *   index of the sorting order from the `column.sorting` initialisation array.
	 *//**
	 * Set the ordering for the table.
	 *
	 * @param {integer} order Column index to sort upon.
	 * @param {string} direction Direction of the sort to be applied (`asc` or `desc`)
	 * @returns {DataTables.Api} this
	 *//**
	 * Set the ordering for the table.
	 *
	 * @param {array} order 1D array of sorting information to be applied.
	 * @param {array} [...] Optional additional sorting conditions
	 * @returns {DataTables.Api} this
	 *//**
	 * Set the ordering for the table.
	 *
	 * @param {array} order 2D array of sorting information to be applied.
	 * @returns {DataTables.Api} this
	 */
	_api_register( 'order()', function ( order, dir ) {
		var ctx = this.context;
	
		if ( order === undefined ) {
			// get
			return ctx.length !== 0 ?
				ctx[0].aaSorting :
				undefined;
		}
	
		// set
		if ( typeof order === 'number' ) {
			// Simple column / direction passed in
			order = [ [ order, dir ] ];
		}
		else if ( order.length && ! Array.isArray( order[0] ) ) {
			// Arguments passed in (list of 1D arrays)
			order = Array.prototype.slice.call( arguments );
		}
		// otherwise a 2D array was passed in
	
		return this.iterator( 'table', function ( settings ) {
			settings.aaSorting = order.slice();
		} );
	} );
	
	
	/**
	 * Attach a sort listener to an element for a given column
	 *
	 * @param {node|jQuery|string} node Identifier for the element(s) to attach the
	 *   listener to. This can take the form of a single DOM node, a jQuery
	 *   collection of nodes or a jQuery selector which will identify the node(s).
	 * @param {integer} column the column that a click on this node will sort on
	 * @param {function} [callback] callback function when sort is run
	 * @returns {DataTables.Api} this
	 */
	_api_register( 'order.listener()', function ( node, column, callback ) {
		return this.iterator( 'table', function ( settings ) {
			_fnSortAttachListener( settings, node, column, callback );
		} );
	} );
	
	
	_api_register( 'order.fixed()', function ( set ) {
		if ( ! set ) {
			var ctx = this.context;
			var fixed = ctx.length ?
				ctx[0].aaSortingFixed :
				undefined;
	
			return Array.isArray( fixed ) ?
				{ pre: fixed } :
				fixed;
		}
	
		return this.iterator( 'table', function ( settings ) {
			settings.aaSortingFixed = $.extend( true, {}, set );
		} );
	} );
	
	
	// Order by the selected column(s)
	_api_register( [
		'columns().order()',
		'column().order()'
	], function ( dir ) {
		var that = this;
	
		return this.iterator( 'table', function ( settings, i ) {
			var sort = [];
	
			$.each( that[i], function (j, col) {
				sort.push( [ col, dir ] );
			} );
	
			settings.aaSorting = sort;
		} );
	} );
	
	
	
	_api_register( 'search()', function ( input, regex, smart, caseInsen ) {
		var ctx = this.context;
	
		if ( input === undefined ) {
			// get
			return ctx.length !== 0 ?
				ctx[0].oPreviousSearch.sSearch :
				undefined;
		}
	
		// set
		return this.iterator( 'table', function ( settings ) {
			if ( ! settings.oFeatures.bFilter ) {
				return;
			}
	
			_fnFilterComplete( settings, $.extend( {}, settings.oPreviousSearch, {
				"sSearch": input+"",
				"bRegex":  regex === null ? false : regex,
				"bSmart":  smart === null ? true  : smart,
				"bCaseInsensitive": caseInsen === null ? true : caseInsen
			} ), 1 );
		} );
	} );
	
	
	_api_registerPlural(
		'columns().search()',
		'column().search()',
		function ( input, regex, smart, caseInsen ) {
			return this.iterator( 'column', function ( settings, column ) {
				var preSearch = settings.aoPreSearchCols;
	
				if ( input === undefined ) {
					// get
					return preSearch[ column ].sSearch;
				}
	
				// set
				if ( ! settings.oFeatures.bFilter ) {
					return;
				}
	
				$.extend( preSearch[ column ], {
					"sSearch": input+"",
					"bRegex":  regex === null ? false : regex,
					"bSmart":  smart === null ? true  : smart,
					"bCaseInsensitive": caseInsen === null ? true : caseInsen
				} );
	
				_fnFilterComplete( settings, settings.oPreviousSearch, 1 );
			} );
		}
	);
	
	/*
	 * State API methods
	 */
	
	_api_register( 'state()', function () {
		return this.context.length ?
			this.context[0].oSavedState :
			null;
	} );
	
	
	_api_register( 'state.clear()', function () {
		return this.iterator( 'table', function ( settings ) {
			// Save an empty object
			settings.fnStateSaveCallback.call( settings.oInstance, settings, {} );
		} );
	} );
	
	
	_api_register( 'state.loaded()', function () {
		return this.context.length ?
			this.context[0].oLoadedState :
			null;
	} );
	
	
	_api_register( 'state.save()', function () {
		return this.iterator( 'table', function ( settings ) {
			_fnSaveState( settings );
		} );
	} );
	
	
	
	/**
	 * Set the jQuery or window object to be used by DataTables
	 *
	 * @param {*} module Library / container object
	 * @param {string} type Library or container type `lib` or `win`.
	 */
	DataTable.use = function (module, type) {
		if (type === 'lib' || module.fn) {
			$ = module;
		}
		else if (type == 'win' || module.document) {
			window = module;
			document = module.document;
		}
	}
	
	/**
	 * CommonJS factory function pass through. This will check if the arguments
	 * given are a window object or a jQuery object. If so they are set
	 * accordingly.
	 * @param {*} root Window
	 * @param {*} jq jQUery
	 * @returns {boolean} Indicator
	 */
	DataTable.factory = function (root, jq) {
		var is = false;
	
		// Test if the first parameter is a window object
		if (root && root.document) {
			window = root;
			document = root.document;
		}
	
		// Test if the second parameter is a jQuery object
		if (jq && jq.fn && jq.fn.jquery) {
			$ = jq;
			is = true;
		}
	
		return is;
	}
	
	/**
	 * Provide a common method for plug-ins to check the version of DataTables being
	 * used, in order to ensure compatibility.
	 *
	 *  @param {string} version Version string to check for, in the format "X.Y.Z".
	 *    Note that the formats "X" and "X.Y" are also acceptable.
	 *  @returns {boolean} true if this version of DataTables is greater or equal to
	 *    the required version, or false if this version of DataTales is not
	 *    suitable
	 *  @static
	 *  @dtopt API-Static
	 *
	 *  @example
	 *    alert( $.fn.dataTable.versionCheck( '1.9.0' ) );
	 */
	DataTable.versionCheck = DataTable.fnVersionCheck = function( version )
	{
		var aThis = DataTable.version.split('.');
		var aThat = version.split('.');
		var iThis, iThat;
	
		for ( var i=0, iLen=aThat.length ; i<iLen ; i++ ) {
			iThis = parseInt( aThis[i], 10 ) || 0;
			iThat = parseInt( aThat[i], 10 ) || 0;
	
			// Parts are the same, keep comparing
			if (iThis === iThat) {
				continue;
			}
	
			// Parts are different, return immediately
			return iThis > iThat;
		}
	
		return true;
	};
	
	
	/**
	 * Check if a `<table>` node is a DataTable table already or not.
	 *
	 *  @param {node|jquery|string} table Table node, jQuery object or jQuery
	 *      selector for the table to test. Note that if more than more than one
	 *      table is passed on, only the first will be checked
	 *  @returns {boolean} true the table given is a DataTable, or false otherwise
	 *  @static
	 *  @dtopt API-Static
	 *
	 *  @example
	 *    if ( ! $.fn.DataTable.isDataTable( '#example' ) ) {
	 *      $('#example').dataTable();
	 *    }
	 */
	DataTable.isDataTable = DataTable.fnIsDataTable = function ( table )
	{
		var t = $(table).get(0);
		var is = false;
	
		if ( table instanceof DataTable.Api ) {
			return true;
		}
	
		$.each( DataTable.settings, function (i, o) {
			var head = o.nScrollHead ? $('table', o.nScrollHead)[0] : null;
			var foot = o.nScrollFoot ? $('table', o.nScrollFoot)[0] : null;
	
			if ( o.nTable === t || head === t || foot === t ) {
				is = true;
			}
		} );
	
		return is;
	};
	
	
	/**
	 * Get all DataTable tables that have been initialised - optionally you can
	 * select to get only currently visible tables.
	 *
	 *  @param {boolean} [visible=false] Flag to indicate if you want all (default)
	 *    or visible tables only.
	 *  @returns {array} Array of `table` nodes (not DataTable instances) which are
	 *    DataTables
	 *  @static
	 *  @dtopt API-Static
	 *
	 *  @example
	 *    $.each( $.fn.dataTable.tables(true), function () {
	 *      $(table).DataTable().columns.adjust();
	 *    } );
	 */
	DataTable.tables = DataTable.fnTables = function ( visible )
	{
		var api = false;
	
		if ( $.isPlainObject( visible ) ) {
			api = visible.api;
			visible = visible.visible;
		}
	
		var a = $.map( DataTable.settings, function (o) {
			if ( !visible || (visible && $(o.nTable).is(':visible')) ) {
				return o.nTable;
			}
		} );
	
		return api ?
			new _Api( a ) :
			a;
	};
	
	
	/**
	 * Convert from camel case parameters to Hungarian notation. This is made public
	 * for the extensions to provide the same ability as DataTables core to accept
	 * either the 1.9 style Hungarian notation, or the 1.10+ style camelCase
	 * parameters.
	 *
	 *  @param {object} src The model object which holds all parameters that can be
	 *    mapped.
	 *  @param {object} user The object to convert from camel case to Hungarian.
	 *  @param {boolean} force When set to `true`, properties which already have a
	 *    Hungarian value in the `user` object will be overwritten. Otherwise they
	 *    won't be.
	 */
	DataTable.camelToHungarian = _fnCamelToHungarian;
	
	
	
	/**
	 *
	 */
	_api_register( '$()', function ( selector, opts ) {
		var
			rows   = this.rows( opts ).nodes(), // Get all rows
			jqRows = $(rows);
	
		return $( [].concat(
			jqRows.filter( selector ).toArray(),
			jqRows.find( selector ).toArray()
		) );
	} );
	
	
	// jQuery functions to operate on the tables
	$.each( [ 'on', 'one', 'off' ], function (i, key) {
		_api_register( key+'()', function ( /* event, handler */ ) {
			var args = Array.prototype.slice.call(arguments);
	
			// Add the `dt` namespace automatically if it isn't already present
			args[0] = $.map( args[0].split( /\s/ ), function ( e ) {
				return ! e.match(/\.dt\b/) ?
					e+'.dt' :
					e;
				} ).join( ' ' );
	
			var inst = $( this.tables().nodes() );
			inst[key].apply( inst, args );
			return this;
		} );
	} );
	
	
	_api_register( 'clear()', function () {
		return this.iterator( 'table', function ( settings ) {
			_fnClearTable( settings );
		} );
	} );
	
	
	_api_register( 'settings()', function () {
		return new _Api( this.context, this.context );
	} );
	
	
	_api_register( 'init()', function () {
		var ctx = this.context;
		return ctx.length ? ctx[0].oInit : null;
	} );
	
	
	_api_register( 'data()', function () {
		return this.iterator( 'table', function ( settings ) {
			return _pluck( settings.aoData, '_aData' );
		} ).flatten();
	} );
	
	
	_api_register( 'destroy()', function ( remove ) {
		remove = remove || false;
	
		return this.iterator( 'table', function ( settings ) {
			var classes   = settings.oClasses;
			var table     = settings.nTable;
			var tbody     = settings.nTBody;
			var thead     = settings.nTHead;
			var tfoot     = settings.nTFoot;
			var jqTable   = $(table);
			var jqTbody   = $(tbody);
			var jqWrapper = $(settings.nTableWrapper);
			var rows      = $.map( settings.aoData, function (r) { return r.nTr; } );
			var i, ien;
	
			// Flag to note that the table is currently being destroyed - no action
			// should be taken
			settings.bDestroying = true;
	
			// Fire off the destroy callbacks for plug-ins etc
			_fnCallbackFire( settings, "aoDestroyCallback", "destroy", [settings] );
	
			// If not being removed from the document, make all columns visible
			if ( ! remove ) {
				new _Api( settings ).columns().visible( true );
			}
	
			// Blitz all `DT` namespaced events (these are internal events, the
			// lowercase, `dt` events are user subscribed and they are responsible
			// for removing them
			jqWrapper.off('.DT').find(':not(tbody *)').off('.DT');
			$(window).off('.DT-'+settings.sInstance);
	
			// When scrolling we had to break the table up - restore it
			if ( table != thead.parentNode ) {
				jqTable.children('thead').detach();
				jqTable.append( thead );
			}
	
			if ( tfoot && table != tfoot.parentNode ) {
				jqTable.children('tfoot').detach();
				jqTable.append( tfoot );
			}
	
			settings.aaSorting = [];
			settings.aaSortingFixed = [];
			_fnSortingClasses( settings );
	
			$( rows ).removeClass( settings.asStripeClasses.join(' ') );
	
			$('th, td', thead).removeClass( classes.sSortable+' '+
				classes.sSortableAsc+' '+classes.sSortableDesc+' '+classes.sSortableNone
			);
	
			// Add the TR elements back into the table in their original order
			jqTbody.children().detach();
			jqTbody.append( rows );
	
			var orig = settings.nTableWrapper.parentNode;
	
			// Remove the DataTables generated nodes, events and classes
			var removedMethod = remove ? 'remove' : 'detach';
			jqTable[ removedMethod ]();
			jqWrapper[ removedMethod ]();
	
			// If we need to reattach the table to the document
			if ( ! remove && orig ) {
				// insertBefore acts like appendChild if !arg[1]
				orig.insertBefore( table, settings.nTableReinsertBefore );
	
				// Restore the width of the original table - was read from the style property,
				// so we can restore directly to that
				jqTable
					.css( 'width', settings.sDestroyWidth )
					.removeClass( classes.sTable );
	
				// If the were originally stripe classes - then we add them back here.
				// Note this is not fool proof (for example if not all rows had stripe
				// classes - but it's a good effort without getting carried away
				ien = settings.asDestroyStripes.length;
	
				if ( ien ) {
					jqTbody.children().each( function (i) {
						$(this).addClass( settings.asDestroyStripes[i % ien] );
					} );
				}
			}
	
			/* Remove the settings object from the settings array */
			var idx = $.inArray( settings, DataTable.settings );
			if ( idx !== -1 ) {
				DataTable.settings.splice( idx, 1 );
			}
		} );
	} );
	
	
	// Add the `every()` method for rows, columns and cells in a compact form
	$.each( [ 'column', 'row', 'cell' ], function ( i, type ) {
		_api_register( type+'s().every()', function ( fn ) {
			var opts = this.selector.opts;
			var api = this;
	
			return this.iterator( type, function ( settings, arg1, arg2, arg3, arg4 ) {
				// Rows and columns:
				//  arg1 - index
				//  arg2 - table counter
				//  arg3 - loop counter
				//  arg4 - undefined
				// Cells:
				//  arg1 - row index
				//  arg2 - column index
				//  arg3 - table counter
				//  arg4 - loop counter
				fn.call(
					api[ type ](
						arg1,
						type==='cell' ? arg2 : opts,
						type==='cell' ? opts : undefined
					),
					arg1, arg2, arg3, arg4
				);
			} );
		} );
	} );
	
	
	// i18n method for extensions to be able to use the language object from the
	// DataTable
	_api_register( 'i18n()', function ( token, def, plural ) {
		var ctx = this.context[0];
		var resolved = _fnGetObjectDataFn( token )( ctx.oLanguage );
	
		if ( resolved === undefined ) {
			resolved = def;
		}
	
		if ( plural !== undefined && $.isPlainObject( resolved ) ) {
			resolved = resolved[ plural ] !== undefined ?
				resolved[ plural ] :
				resolved._;
		}
	
		return resolved.replace( '%d', plural ); // nb: plural might be undefined,
	} );	
	/**
	 * Version string for plug-ins to check compatibility. Allowed format is
	 * `a.b.c-d` where: a:int, b:int, c:int, d:string(dev|beta|alpha). `d` is used
	 * only for non-release builds. See http://semver.org/ for more information.
	 *  @member
	 *  @type string
	 *  @default Version number
	 */
	DataTable.version = "1.13.4";
	
	/**
	 * Private data store, containing all of the settings objects that are
	 * created for the tables on a given page.
	 *
	 * Note that the `DataTable.settings` object is aliased to
	 * `jQuery.fn.dataTableExt` through which it may be accessed and
	 * manipulated, or `jQuery.fn.dataTable.settings`.
	 *  @member
	 *  @type array
	 *  @default []
	 *  @private
	 */
	DataTable.settings = [];
	
	/**
	 * Object models container, for the various models that DataTables has
	 * available to it. These models define the objects that are used to hold
	 * the active state and configuration of the table.
	 *  @namespace
	 */
	DataTable.models = {};
	
	
	
	/**
	 * Template object for the way in which DataTables holds information about
	 * search information for the global filter and individual column filters.
	 *  @namespace
	 */
	DataTable.models.oSearch = {
		/**
		 * Flag to indicate if the filtering should be case insensitive or not
		 *  @type boolean
		 *  @default true
		 */
		"bCaseInsensitive": true,
	
		/**
		 * Applied search term
		 *  @type string
		 *  @default <i>Empty string</i>
		 */
		"sSearch": "",
	
		/**
		 * Flag to indicate if the search term should be interpreted as a
		 * regular expression (true) or not (false) and therefore and special
		 * regex characters escaped.
		 *  @type boolean
		 *  @default false
		 */
		"bRegex": false,
	
		/**
		 * Flag to indicate if DataTables is to use its smart filtering or not.
		 *  @type boolean
		 *  @default true
		 */
		"bSmart": true,
	
		/**
		 * Flag to indicate if DataTables should only trigger a search when
		 * the return key is pressed.
		 *  @type boolean
		 *  @default false
		 */
		"return": false
	};
	
	
	
	
	/**
	 * Template object for the way in which DataTables holds information about
	 * each individual row. This is the object format used for the settings
	 * aoData array.
	 *  @namespace
	 */
	DataTable.models.oRow = {
		/**
		 * TR element for the row
		 *  @type node
		 *  @default null
		 */
		"nTr": null,
	
		/**
		 * Array of TD elements for each row. This is null until the row has been
		 * created.
		 *  @type array nodes
		 *  @default []
		 */
		"anCells": null,
	
		/**
		 * Data object from the original data source for the row. This is either
		 * an array if using the traditional form of DataTables, or an object if
		 * using mData options. The exact type will depend on the passed in
		 * data from the data source, or will be an array if using DOM a data
		 * source.
		 *  @type array|object
		 *  @default []
		 */
		"_aData": [],
	
		/**
		 * Sorting data cache - this array is ostensibly the same length as the
		 * number of columns (although each index is generated only as it is
		 * needed), and holds the data that is used for sorting each column in the
		 * row. We do this cache generation at the start of the sort in order that
		 * the formatting of the sort data need be done only once for each cell
		 * per sort. This array should not be read from or written to by anything
		 * other than the master sorting methods.
		 *  @type array
		 *  @default null
		 *  @private
		 */
		"_aSortData": null,
	
		/**
		 * Per cell filtering data cache. As per the sort data cache, used to
		 * increase the performance of the filtering in DataTables
		 *  @type array
		 *  @default null
		 *  @private
		 */
		"_aFilterData": null,
	
		/**
		 * Filtering data cache. This is the same as the cell filtering cache, but
		 * in this case a string rather than an array. This is easily computed with
		 * a join on `_aFilterData`, but is provided as a cache so the join isn't
		 * needed on every search (memory traded for performance)
		 *  @type array
		 *  @default null
		 *  @private
		 */
		"_sFilterRow": null,
	
		/**
		 * Cache of the class name that DataTables has applied to the row, so we
		 * can quickly look at this variable rather than needing to do a DOM check
		 * on className for the nTr property.
		 *  @type string
		 *  @default <i>Empty string</i>
		 *  @private
		 */
		"_sRowStripe": "",
	
		/**
		 * Denote if the original data source was from the DOM, or the data source
		 * object. This is used for invalidating data, so DataTables can
		 * automatically read data from the original source, unless uninstructed
		 * otherwise.
		 *  @type string
		 *  @default null
		 *  @private
		 */
		"src": null,
	
		/**
		 * Index in the aoData array. This saves an indexOf lookup when we have the
		 * object, but want to know the index
		 *  @type integer
		 *  @default -1
		 *  @private
		 */
		"idx": -1
	};
	
	
	/**
	 * Template object for the column information object in DataTables. This object
	 * is held in the settings aoColumns array and contains all the information that
	 * DataTables needs about each individual column.
	 *
	 * Note that this object is related to {@link DataTable.defaults.column}
	 * but this one is the internal data store for DataTables's cache of columns.
	 * It should NOT be manipulated outside of DataTables. Any configuration should
	 * be done through the initialisation options.
	 *  @namespace
	 */
	DataTable.models.oColumn = {
		/**
		 * Column index. This could be worked out on-the-fly with $.inArray, but it
		 * is faster to just hold it as a variable
		 *  @type integer
		 *  @default null
		 */
		"idx": null,
	
		/**
		 * A list of the columns that sorting should occur on when this column
		 * is sorted. That this property is an array allows multi-column sorting
		 * to be defined for a column (for example first name / last name columns
		 * would benefit from this). The values are integers pointing to the
		 * columns to be sorted on (typically it will be a single integer pointing
		 * at itself, but that doesn't need to be the case).
		 *  @type array
		 */
		"aDataSort": null,
	
		/**
		 * Define the sorting directions that are applied to the column, in sequence
		 * as the column is repeatedly sorted upon - i.e. the first value is used
		 * as the sorting direction when the column if first sorted (clicked on).
		 * Sort it again (click again) and it will move on to the next index.
		 * Repeat until loop.
		 *  @type array
		 */
		"asSorting": null,
	
		/**
		 * Flag to indicate if the column is searchable, and thus should be included
		 * in the filtering or not.
		 *  @type boolean
		 */
		"bSearchable": null,
	
		/**
		 * Flag to indicate if the column is sortable or not.
		 *  @type boolean
		 */
		"bSortable": null,
	
		/**
		 * Flag to indicate if the column is currently visible in the table or not
		 *  @type boolean
		 */
		"bVisible": null,
	
		/**
		 * Store for manual type assignment using the `column.type` option. This
		 * is held in store so we can manipulate the column's `sType` property.
		 *  @type string
		 *  @default null
		 *  @private
		 */
		"_sManualType": null,
	
		/**
		 * Flag to indicate if HTML5 data attributes should be used as the data
		 * source for filtering or sorting. True is either are.
		 *  @type boolean
		 *  @default false
		 *  @private
		 */
		"_bAttrSrc": false,
	
		/**
		 * Developer definable function that is called whenever a cell is created (Ajax source,
		 * etc) or processed for input (DOM source). This can be used as a compliment to mRender
		 * allowing you to modify the DOM element (add background colour for example) when the
		 * element is available.
		 *  @type function
		 *  @param {element} nTd The TD node that has been created
		 *  @param {*} sData The Data for the cell
		 *  @param {array|object} oData The data for the whole row
		 *  @param {int} iRow The row index for the aoData data store
		 *  @default null
		 */
		"fnCreatedCell": null,
	
		/**
		 * Function to get data from a cell in a column. You should <b>never</b>
		 * access data directly through _aData internally in DataTables - always use
		 * the method attached to this property. It allows mData to function as
		 * required. This function is automatically assigned by the column
		 * initialisation method
		 *  @type function
		 *  @param {array|object} oData The data array/object for the array
		 *    (i.e. aoData[]._aData)
		 *  @param {string} sSpecific The specific data type you want to get -
		 *    'display', 'type' 'filter' 'sort'
		 *  @returns {*} The data for the cell from the given row's data
		 *  @default null
		 */
		"fnGetData": null,
	
		/**
		 * Function to set data for a cell in the column. You should <b>never</b>
		 * set the data directly to _aData internally in DataTables - always use
		 * this method. It allows mData to function as required. This function
		 * is automatically assigned by the column initialisation method
		 *  @type function
		 *  @param {array|object} oData The data array/object for the array
		 *    (i.e. aoData[]._aData)
		 *  @param {*} sValue Value to set
		 *  @default null
		 */
		"fnSetData": null,
	
		/**
		 * Property to read the value for the cells in the column from the data
		 * source array / object. If null, then the default content is used, if a
		 * function is given then the return from the function is used.
		 *  @type function|int|string|null
		 *  @default null
		 */
		"mData": null,
	
		/**
		 * Partner property to mData which is used (only when defined) to get
		 * the data - i.e. it is basically the same as mData, but without the
		 * 'set' option, and also the data fed to it is the result from mData.
		 * This is the rendering method to match the data method of mData.
		 *  @type function|int|string|null
		 *  @default null
		 */
		"mRender": null,
	
		/**
		 * Unique header TH/TD element for this column - this is what the sorting
		 * listener is attached to (if sorting is enabled.)
		 *  @type node
		 *  @default null
		 */
		"nTh": null,
	
		/**
		 * Unique footer TH/TD element for this column (if there is one). Not used
		 * in DataTables as such, but can be used for plug-ins to reference the
		 * footer for each column.
		 *  @type node
		 *  @default null
		 */
		"nTf": null,
	
		/**
		 * The class to apply to all TD elements in the table's TBODY for the column
		 *  @type string
		 *  @default null
		 */
		"sClass": null,
	
		/**
		 * When DataTables calculates the column widths to assign to each column,
		 * it finds the longest string in each column and then constructs a
		 * temporary table and reads the widths from that. The problem with this
		 * is that "mmm" is much wider then "iiii", but the latter is a longer
		 * string - thus the calculation can go wrong (doing it properly and putting
		 * it into an DOM object and measuring that is horribly(!) slow). Thus as
		 * a "work around" we provide this option. It will append its value to the
		 * text that is found to be the longest string for the column - i.e. padding.
		 *  @type string
		 */
		"sContentPadding": null,
	
		/**
		 * Allows a default value to be given for a column's data, and will be used
		 * whenever a null data source is encountered (this can be because mData
		 * is set to null, or because the data source itself is null).
		 *  @type string
		 *  @default null
		 */
		"sDefaultContent": null,
	
		/**
		 * Name for the column, allowing reference to the column by name as well as
		 * by index (needs a lookup to work by name).
		 *  @type string
		 */
		"sName": null,
	
		/**
		 * Custom sorting data type - defines which of the available plug-ins in
		 * afnSortData the custom sorting will use - if any is defined.
		 *  @type string
		 *  @default std
		 */
		"sSortDataType": 'std',
	
		/**
		 * Class to be applied to the header element when sorting on this column
		 *  @type string
		 *  @default null
		 */
		"sSortingClass": null,
	
		/**
		 * Class to be applied to the header element when sorting on this column -
		 * when jQuery UI theming is used.
		 *  @type string
		 *  @default null
		 */
		"sSortingClassJUI": null,
	
		/**
		 * Title of the column - what is seen in the TH element (nTh).
		 *  @type string
		 */
		"sTitle": null,
	
		/**
		 * Column sorting and filtering type
		 *  @type string
		 *  @default null
		 */
		"sType": null,
	
		/**
		 * Width of the column
		 *  @type string
		 *  @default null
		 */
		"sWidth": null,
	
		/**
		 * Width of the column when it was first "encountered"
		 *  @type string
		 *  @default null
		 */
		"sWidthOrig": null
	};
	
	
	/*
	 * Developer note: The properties of the object below are given in Hungarian
	 * notation, that was used as the interface for DataTables prior to v1.10, however
	 * from v1.10 onwards the primary interface is camel case. In order to avoid
	 * breaking backwards compatibility utterly with this change, the Hungarian
	 * version is still, internally the primary interface, but is is not documented
	 * - hence the @name tags in each doc comment. This allows a Javascript function
	 * to create a map from Hungarian notation to camel case (going the other direction
	 * would require each property to be listed, which would add around 3K to the size
	 * of DataTables, while this method is about a 0.5K hit).
	 *
	 * Ultimately this does pave the way for Hungarian notation to be dropped
	 * completely, but that is a massive amount of work and will break current
	 * installs (therefore is on-hold until v2).
	 */
	
	/**
	 * Initialisation options that can be given to DataTables at initialisation
	 * time.
	 *  @namespace
	 */
	DataTable.defaults = {
		/**
		 * An array of data to use for the table, passed in at initialisation which
		 * will be used in preference to any data which is already in the DOM. This is
		 * particularly useful for constructing tables purely in Javascript, for
		 * example with a custom Ajax call.
		 *  @type array
		 *  @default null
		 *
		 *  @dtopt Option
		 *  @name DataTable.defaults.data
		 *
		 *  @example
		 *    // Using a 2D array data source
		 *    $(document).ready( function () {
		 *      $('#example').dataTable( {
		 *        "data": [
		 *          ['Trident', 'Internet Explorer 4.0', 'Win 95+', 4, 'X'],
		 *          ['Trident', 'Internet Explorer 5.0', 'Win 95+', 5, 'C'],
		 *        ],
		 *        "columns": [
		 *          { "title": "Engine" },
		 *          { "title": "Browser" },
		 *          { "title": "Platform" },
		 *          { "title": "Version" },
		 *          { "title": "Grade" }
		 *        ]
		 *      } );
		 *    } );
		 *
		 *  @example
		 *    // Using an array of objects as a data source (`data`)
		 *    $(document).ready( function () {
		 *      $('#example').dataTable( {
		 *        "data": [
		 *          {
		 *            "engine":   "Trident",
		 *            "browser":  "Internet Explorer 4.0",
		 *            "platform": "Win 95+",
		 *            "version":  4,
		 *            "grade":    "X"
		 *          },
		 *          {
		 *            "engine":   "Trident",
		 *            "browser":  "Internet Explorer 5.0",
		 *            "platform": "Win 95+",
		 *            "version":  5,
		 *            "grade":    "C"
		 *          }
		 *        ],
		 *        "columns": [
		 *          { "title": "Engine",   "data": "engine" },
		 *          { "title": "Browser",  "data": "browser" },
		 *          { "title": "Platform", "data": "platform" },
		 *          { "title": "Version",  "data": "version" },
		 *          { "title": "Grade",    "data": "grade" }
		 *        ]
		 *      } );
		 *    } );
		 */
		"aaData": null,
	
	
		/**
		 * If ordering is enabled, then DataTables will perform a first pass sort on
		 * initialisation. You can define which column(s) the sort is performed
		 * upon, and the sorting direction, with this variable. The `sorting` array
		 * should contain an array for each column to be sorted initially containing
		 * the column's index and a direction string ('asc' or 'desc').
		 *  @type array
		 *  @default [[0,'asc']]
		 *
		 *  @dtopt Option
		 *  @name DataTable.defaults.order
		 *
		 *  @example
		 *    // Sort by 3rd column first, and then 4th column
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "order": [[2,'asc'], [3,'desc']]
		 *      } );
		 *    } );
		 *
		 *    // No initial sorting
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "order": []
		 *      } );
		 *    } );
		 */
		"aaSorting": [[0,'asc']],
	
	
		/**
		 * This parameter is basically identical to the `sorting` parameter, but
		 * cannot be overridden by user interaction with the table. What this means
		 * is that you could have a column (visible or hidden) which the sorting
		 * will always be forced on first - any sorting after that (from the user)
		 * will then be performed as required. This can be useful for grouping rows
		 * together.
		 *  @type array
		 *  @default null
		 *
		 *  @dtopt Option
		 *  @name DataTable.defaults.orderFixed
		 *
		 *  @example
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "orderFixed": [[0,'asc']]
		 *      } );
		 *    } )
		 */
		"aaSortingFixed": [],
	
	
		/**
		 * DataTables can be instructed to load data to display in the table from a
		 * Ajax source. This option defines how that Ajax call is made and where to.
		 *
		 * The `ajax` property has three different modes of operation, depending on
		 * how it is defined. These are:
		 *
		 * * `string` - Set the URL from where the data should be loaded from.
		 * * `object` - Define properties for `jQuery.ajax`.
		 * * `function` - Custom data get function
		 *
		 * `string`
		 * --------
		 *
		 * As a string, the `ajax` property simply defines the URL from which
		 * DataTables will load data.
		 *
		 * `object`
		 * --------
		 *
		 * As an object, the parameters in the object are passed to
		 * [jQuery.ajax](http://api.jquery.com/jQuery.ajax/) allowing fine control
		 * of the Ajax request. DataTables has a number of default parameters which
		 * you can override using this option. Please refer to the jQuery
		 * documentation for a full description of the options available, although
		 * the following parameters provide additional options in DataTables or
		 * require special consideration:
		 *
		 * * `data` - As with jQuery, `data` can be provided as an object, but it
		 *   can also be used as a function to manipulate the data DataTables sends
		 *   to the server. The function takes a single parameter, an object of
		 *   parameters with the values that DataTables has readied for sending. An
		 *   object may be returned which will be merged into the DataTables
		 *   defaults, or you can add the items to the object that was passed in and
		 *   not return anything from the function. This supersedes `fnServerParams`
		 *   from DataTables 1.9-.
		 *
		 * * `dataSrc` - By default DataTables will look for the property `data` (or
		 *   `aaData` for compatibility with DataTables 1.9-) when obtaining data
		 *   from an Ajax source or for server-side processing - this parameter
		 *   allows that property to be changed. You can use Javascript dotted
		 *   object notation to get a data source for multiple levels of nesting, or
		 *   it my be used as a function. As a function it takes a single parameter,
		 *   the JSON returned from the server, which can be manipulated as
		 *   required, with the returned value being that used by DataTables as the
		 *   data source for the table. This supersedes `sAjaxDataProp` from
		 *   DataTables 1.9-.
		 *
		 * * `success` - Should not be overridden it is used internally in
		 *   DataTables. To manipulate / transform the data returned by the server
		 *   use `ajax.dataSrc`, or use `ajax` as a function (see below).
		 *
		 * `function`
		 * ----------
		 *
		 * As a function, making the Ajax call is left up to yourself allowing
		 * complete control of the Ajax request. Indeed, if desired, a method other
		 * than Ajax could be used to obtain the required data, such as Web storage
		 * or an AIR database.
		 *
		 * The function is given four parameters and no return is required. The
		 * parameters are:
		 *
		 * 1. _object_ - Data to send to the server
		 * 2. _function_ - Callback function that must be executed when the required
		 *    data has been obtained. That data should be passed into the callback
		 *    as the only parameter
		 * 3. _object_ - DataTables settings object for the table
		 *
		 * Note that this supersedes `fnServerData` from DataTables 1.9-.
		 *
		 *  @type string|object|function
		 *  @default null
		 *
		 *  @dtopt Option
		 *  @name DataTable.defaults.ajax
		 *  @since 1.10.0
		 *
		 * @example
		 *   // Get JSON data from a file via Ajax.
		 *   // Note DataTables expects data in the form `{ data: [ ...data... ] }` by default).
		 *   $('#example').dataTable( {
		 *     "ajax": "data.json"
		 *   } );
		 *
		 * @example
		 *   // Get JSON data from a file via Ajax, using `dataSrc` to change
		 *   // `data` to `tableData` (i.e. `{ tableData: [ ...data... ] }`)
		 *   $('#example').dataTable( {
		 *     "ajax": {
		 *       "url": "data.json",
		 *       "dataSrc": "tableData"
		 *     }
		 *   } );
		 *
		 * @example
		 *   // Get JSON data from a file via Ajax, using `dataSrc` to read data
		 *   // from a plain array rather than an array in an object
		 *   $('#example').dataTable( {
		 *     "ajax": {
		 *       "url": "data.json",
		 *       "dataSrc": ""
		 *     }
		 *   } );
		 *
		 * @example
		 *   // Manipulate the data returned from the server - add a link to data
		 *   // (note this can, should, be done using `render` for the column - this
		 *   // is just a simple example of how the data can be manipulated).
		 *   $('#example').dataTable( {
		 *     "ajax": {
		 *       "url": "data.json",
		 *       "dataSrc": function ( json ) {
		 *         for ( var i=0, ien=json.length ; i<ien ; i++ ) {
		 *           json[i][0] = '<a href="/message/'+json[i][0]+'>View message</a>';
		 *         }
		 *         return json;
		 *       }
		 *     }
		 *   } );
		 *
		 * @example
		 *   // Add data to the request
		 *   $('#example').dataTable( {
		 *     "ajax": {
		 *       "url": "data.json",
		 *       "data": function ( d ) {
		 *         return {
		 *           "extra_search": $('#extra').val()
		 *         };
		 *       }
		 *     }
		 *   } );
		 *
		 * @example
		 *   // Send request as POST
		 *   $('#example').dataTable( {
		 *     "ajax": {
		 *       "url": "data.json",
		 *       "type": "POST"
		 *     }
		 *   } );
		 *
		 * @example
		 *   // Get the data from localStorage (could interface with a form for
		 *   // adding, editing and removing rows).
		 *   $('#example').dataTable( {
		 *     "ajax": function (data, callback, settings) {
		 *       callback(
		 *         JSON.parse( localStorage.getItem('dataTablesData') )
		 *       );
		 *     }
		 *   } );
		 */
		"ajax": null,
	
	
		/**
		 * This parameter allows you to readily specify the entries in the length drop
		 * down menu that DataTables shows when pagination is enabled. It can be
		 * either a 1D array of options which will be used for both the displayed
		 * option and the value, or a 2D array which will use the array in the first
		 * position as the value, and the array in the second position as the
		 * displayed options (useful for language strings such as 'All').
		 *
		 * Note that the `pageLength` property will be automatically set to the
		 * first value given in this array, unless `pageLength` is also provided.
		 *  @type array
		 *  @default [ 10, 25, 50, 100 ]
		 *
		 *  @dtopt Option
		 *  @name DataTable.defaults.lengthMenu
		 *
		 *  @example
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "lengthMenu": [[10, 25, 50, -1], [10, 25, 50, "All"]]
		 *      } );
		 *    } );
		 */
		"aLengthMenu": [ 10, 25, 50, 100 ],
	
	
		/**
		 * The `columns` option in the initialisation parameter allows you to define
		 * details about the way individual columns behave. For a full list of
		 * column options that can be set, please see
		 * {@link DataTable.defaults.column}. Note that if you use `columns` to
		 * define your columns, you must have an entry in the array for every single
		 * column that you have in your table (these can be null if you don't which
		 * to specify any options).
		 *  @member
		 *
		 *  @name DataTable.defaults.column
		 */
		"aoColumns": null,
	
		/**
		 * Very similar to `columns`, `columnDefs` allows you to target a specific
		 * column, multiple columns, or all columns, using the `targets` property of
		 * each object in the array. This allows great flexibility when creating
		 * tables, as the `columnDefs` arrays can be of any length, targeting the
		 * columns you specifically want. `columnDefs` may use any of the column
		 * options available: {@link DataTable.defaults.column}, but it _must_
		 * have `targets` defined in each object in the array. Values in the `targets`
		 * array may be:
		 *   <ul>
		 *     <li>a string - class name will be matched on the TH for the column</li>
		 *     <li>0 or a positive integer - column index counting from the left</li>
		 *     <li>a negative integer - column index counting from the right</li>
		 *     <li>the string "_all" - all columns (i.e. assign a default)</li>
		 *   </ul>
		 *  @member
		 *
		 *  @name DataTable.defaults.columnDefs
		 */
		"aoColumnDefs": null,
	
	
		/**
		 * Basically the same as `search`, this parameter defines the individual column
		 * filtering state at initialisation time. The array must be of the same size
		 * as the number of columns, and each element be an object with the parameters
		 * `search` and `escapeRegex` (the latter is optional). 'null' is also
		 * accepted and the default will be used.
		 *  @type array
		 *  @default []
		 *
		 *  @dtopt Option
		 *  @name DataTable.defaults.searchCols
		 *
		 *  @example
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "searchCols": [
		 *          null,
		 *          { "search": "My filter" },
		 *          null,
		 *          { "search": "^[0-9]", "escapeRegex": false }
		 *        ]
		 *      } );
		 *    } )
		 */
		"aoSearchCols": [],
	
	
		/**
		 * An array of CSS classes that should be applied to displayed rows. This
		 * array may be of any length, and DataTables will apply each class
		 * sequentially, looping when required.
		 *  @type array
		 *  @default null <i>Will take the values determined by the `oClasses.stripe*`
		 *    options</i>
		 *
		 *  @dtopt Option
		 *  @name DataTable.defaults.stripeClasses
		 *
		 *  @example
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "stripeClasses": [ 'strip1', 'strip2', 'strip3' ]
		 *      } );
		 *    } )
		 */
		"asStripeClasses": null,
	
	
		/**
		 * Enable or disable automatic column width calculation. This can be disabled
		 * as an optimisation (it takes some time to calculate the widths) if the
		 * tables widths are passed in using `columns`.
		 *  @type boolean
		 *  @default true
		 *
		 *  @dtopt Features
		 *  @name DataTable.defaults.autoWidth
		 *
		 *  @example
		 *    $(document).ready( function () {
		 *      $('#example').dataTable( {
		 *        "autoWidth": false
		 *      } );
		 *    } );
		 */
		"bAutoWidth": true,
	
	
		/**
		 * Deferred rendering can provide DataTables with a huge speed boost when you
		 * are using an Ajax or JS data source for the table. This option, when set to
		 * true, will cause DataTables to defer the creation of the table elements for
		 * each row until they are needed for a draw - saving a significant amount of
		 * time.
		 *  @type boolean
		 *  @default false
		 *
		 *  @dtopt Features
		 *  @name DataTable.defaults.deferRender
		 *
		 *  @example
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "ajax": "sources/arrays.txt",
		 *        "deferRender": true
		 *      } );
		 *    } );
		 */
		"bDeferRender": false,
	
	
		/**
		 * Replace a DataTable which matches the given selector and replace it with
		 * one which has the properties of the new initialisation object passed. If no
		 * table matches the selector, then the new DataTable will be constructed as
		 * per normal.
		 *  @type boolean
		 *  @default false
		 *
		 *  @dtopt Options
		 *  @name DataTable.defaults.destroy
		 *
		 *  @example
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "srollY": "200px",
		 *        "paginate": false
		 *      } );
		 *
		 *      // Some time later....
		 *      $('#example').dataTable( {
		 *        "filter": false,
		 *        "destroy": true
		 *      } );
		 *    } );
		 */
		"bDestroy": false,
	
	
		/**
		 * Enable or disable filtering of data. Filtering in DataTables is "smart" in
		 * that it allows the end user to input multiple words (space separated) and
		 * will match a row containing those words, even if not in the order that was
		 * specified (this allow matching across multiple columns). Note that if you
		 * wish to use filtering in DataTables this must remain 'true' - to remove the
		 * default filtering input box and retain filtering abilities, please use
		 * {@link DataTable.defaults.dom}.
		 *  @type boolean
		 *  @default true
		 *
		 *  @dtopt Features
		 *  @name DataTable.defaults.searching
		 *
		 *  @example
		 *    $(document).ready( function () {
		 *      $('#example').dataTable( {
		 *        "searching": false
		 *      } );
		 *    } );
		 */
		"bFilter": true,
	
	
		/**
		 * Enable or disable the table information display. This shows information
		 * about the data that is currently visible on the page, including information
		 * about filtered data if that action is being performed.
		 *  @type boolean
		 *  @default true
		 *
		 *  @dtopt Features
		 *  @name DataTable.defaults.info
		 *
		 *  @example
		 *    $(document).ready( function () {
		 *      $('#example').dataTable( {
		 *        "info": false
		 *      } );
		 *    } );
		 */
		"bInfo": true,
	
	
		/**
		 * Allows the end user to select the size of a formatted page from a select
		 * menu (sizes are 10, 25, 50 and 100). Requires pagination (`paginate`).
		 *  @type boolean
		 *  @default true
		 *
		 *  @dtopt Features
		 *  @name DataTable.defaults.lengthChange
		 *
		 *  @example
		 *    $(document).ready( function () {
		 *      $('#example').dataTable( {
		 *        "lengthChange": false
		 *      } );
		 *    } );
		 */
		"bLengthChange": true,
	
	
		/**
		 * Enable or disable pagination.
		 *  @type boolean
		 *  @default true
		 *
		 *  @dtopt Features
		 *  @name DataTable.defaults.paging
		 *
		 *  @example
		 *    $(document).ready( function () {
		 *      $('#example').dataTable( {
		 *        "paging": false
		 *      } );
		 *    } );
		 */
		"bPaginate": true,
	
	
		/**
		 * Enable or disable the display of a 'processing' indicator when the table is
		 * being processed (e.g. a sort). This is particularly useful for tables with
		 * large amounts of data where it can take a noticeable amount of time to sort
		 * the entries.
		 *  @type boolean
		 *  @default false
		 *
		 *  @dtopt Features
		 *  @name DataTable.defaults.processing
		 *
		 *  @example
		 *    $(document).ready( function () {
		 *      $('#example').dataTable( {
		 *        "processing": true
		 *      } );
		 *    } );
		 */
		"bProcessing": false,
	
	
		/**
		 * Retrieve the DataTables object for the given selector. Note that if the
		 * table has already been initialised, this parameter will cause DataTables
		 * to simply return the object that has already been set up - it will not take
		 * account of any changes you might have made to the initialisation object
		 * passed to DataTables (setting this parameter to true is an acknowledgement
		 * that you understand this). `destroy` can be used to reinitialise a table if
		 * you need.
		 *  @type boolean
		 *  @default false
		 *
		 *  @dtopt Options
		 *  @name DataTable.defaults.retrieve
		 *
		 *  @example
		 *    $(document).ready( function() {
		 *      initTable();
		 *      tableActions();
		 *    } );
		 *
		 *    function initTable ()
		 *    {
		 *      return $('#example').dataTable( {
		 *        "scrollY": "200px",
		 *        "paginate": false,
		 *        "retrieve": true
		 *      } );
		 *    }
		 *
		 *    function tableActions ()
		 *    {
		 *      var table = initTable();
		 *      // perform API operations with oTable
		 *    }
		 */
		"bRetrieve": false,
	
	
		/**
		 * When vertical (y) scrolling is enabled, DataTables will force the height of
		 * the table's viewport to the given height at all times (useful for layout).
		 * However, this can look odd when filtering data down to a small data set,
		 * and the footer is left "floating" further down. This parameter (when
		 * enabled) will cause DataTables to collapse the table's viewport down when
		 * the result set will fit within the given Y height.
		 *  @type boolean
		 *  @default false
		 *
		 *  @dtopt Options
		 *  @name DataTable.defaults.scrollCollapse
		 *
		 *  @example
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "scrollY": "200",
		 *        "scrollCollapse": true
		 *      } );
		 *    } );
		 */
		"bScrollCollapse": false,
	
	
		/**
		 * Configure DataTables to use server-side processing. Note that the
		 * `ajax` parameter must also be given in order to give DataTables a
		 * source to obtain the required data for each draw.
		 *  @type boolean
		 *  @default false
		 *
		 *  @dtopt Features
		 *  @dtopt Server-side
		 *  @name DataTable.defaults.serverSide
		 *
		 *  @example
		 *    $(document).ready( function () {
		 *      $('#example').dataTable( {
		 *        "serverSide": true,
		 *        "ajax": "xhr.php"
		 *      } );
		 *    } );
		 */
		"bServerSide": false,
	
	
		/**
		 * Enable or disable sorting of columns. Sorting of individual columns can be
		 * disabled by the `sortable` option for each column.
		 *  @type boolean
		 *  @default true
		 *
		 *  @dtopt Features
		 *  @name DataTable.defaults.ordering
		 *
		 *  @example
		 *    $(document).ready( function () {
		 *      $('#example').dataTable( {
		 *        "ordering": false
		 *      } );
		 *    } );
		 */
		"bSort": true,
	
	
		/**
		 * Enable or display DataTables' ability to sort multiple columns at the
		 * same time (activated by shift-click by the user).
		 *  @type boolean
		 *  @default true
		 *
		 *  @dtopt Options
		 *  @name DataTable.defaults.orderMulti
		 *
		 *  @example
		 *    // Disable multiple column sorting ability
		 *    $(document).ready( function () {
		 *      $('#example').dataTable( {
		 *        "orderMulti": false
		 *      } );
		 *    } );
		 */
		"bSortMulti": true,
	
	
		/**
		 * Allows control over whether DataTables should use the top (true) unique
		 * cell that is found for a single column, or the bottom (false - default).
		 * This is useful when using complex headers.
		 *  @type boolean
		 *  @default false
		 *
		 *  @dtopt Options
		 *  @name DataTable.defaults.orderCellsTop
		 *
		 *  @example
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "orderCellsTop": true
		 *      } );
		 *    } );
		 */
		"bSortCellsTop": false,
	
	
		/**
		 * Enable or disable the addition of the classes `sorting\_1`, `sorting\_2` and
		 * `sorting\_3` to the columns which are currently being sorted on. This is
		 * presented as a feature switch as it can increase processing time (while
		 * classes are removed and added) so for large data sets you might want to
		 * turn this off.
		 *  @type boolean
		 *  @default true
		 *
		 *  @dtopt Features
		 *  @name DataTable.defaults.orderClasses
		 *
		 *  @example
		 *    $(document).ready( function () {
		 *      $('#example').dataTable( {
		 *        "orderClasses": false
		 *      } );
		 *    } );
		 */
		"bSortClasses": true,
	
	
		/**
		 * Enable or disable state saving. When enabled HTML5 `localStorage` will be
		 * used to save table display information such as pagination information,
		 * display length, filtering and sorting. As such when the end user reloads
		 * the page the display display will match what thy had previously set up.
		 *
		 * Due to the use of `localStorage` the default state saving is not supported
		 * in IE6 or 7. If state saving is required in those browsers, use
		 * `stateSaveCallback` to provide a storage solution such as cookies.
		 *  @type boolean
		 *  @default false
		 *
		 *  @dtopt Features
		 *  @name DataTable.defaults.stateSave
		 *
		 *  @example
		 *    $(document).ready( function () {
		 *      $('#example').dataTable( {
		 *        "stateSave": true
		 *      } );
		 *    } );
		 */
		"bStateSave": false,
	
	
		/**
		 * This function is called when a TR element is created (and all TD child
		 * elements have been inserted), or registered if using a DOM source, allowing
		 * manipulation of the TR element (adding classes etc).
		 *  @type function
		 *  @param {node} row "TR" element for the current row
		 *  @param {array} data Raw data array for this row
		 *  @param {int} dataIndex The index of this row in the internal aoData array
		 *
		 *  @dtopt Callbacks
		 *  @name DataTable.defaults.createdRow
		 *
		 *  @example
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "createdRow": function( row, data, dataIndex ) {
		 *          // Bold the grade for all 'A' grade browsers
		 *          if ( data[4] == "A" )
		 *          {
		 *            $('td:eq(4)', row).html( '<b>A</b>' );
		 *          }
		 *        }
		 *      } );
		 *    } );
		 */
		"fnCreatedRow": null,
	
	
		/**
		 * This function is called on every 'draw' event, and allows you to
		 * dynamically modify any aspect you want about the created DOM.
		 *  @type function
		 *  @param {object} settings DataTables settings object
		 *
		 *  @dtopt Callbacks
		 *  @name DataTable.defaults.drawCallback
		 *
		 *  @example
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "drawCallback": function( settings ) {
		 *          alert( 'DataTables has redrawn the table' );
		 *        }
		 *      } );
		 *    } );
		 */
		"fnDrawCallback": null,
	
	
		/**
		 * Identical to fnHeaderCallback() but for the table footer this function
		 * allows you to modify the table footer on every 'draw' event.
		 *  @type function
		 *  @param {node} foot "TR" element for the footer
		 *  @param {array} data Full table data (as derived from the original HTML)
		 *  @param {int} start Index for the current display starting point in the
		 *    display array
		 *  @param {int} end Index for the current display ending point in the
		 *    display array
		 *  @param {array int} display Index array to translate the visual position
		 *    to the full data array
		 *
		 *  @dtopt Callbacks
		 *  @name DataTable.defaults.footerCallback
		 *
		 *  @example
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "footerCallback": function( tfoot, data, start, end, display ) {
		 *          tfoot.getElementsByTagName('th')[0].innerHTML = "Starting index is "+start;
		 *        }
		 *      } );
		 *    } )
		 */
		"fnFooterCallback": null,
	
	
		/**
		 * When rendering large numbers in the information element for the table
		 * (i.e. "Showing 1 to 10 of 57 entries") DataTables will render large numbers
		 * to have a comma separator for the 'thousands' units (e.g. 1 million is
		 * rendered as "1,000,000") to help readability for the end user. This
		 * function will override the default method DataTables uses.
		 *  @type function
		 *  @member
		 *  @param {int} toFormat number to be formatted
		 *  @returns {string} formatted string for DataTables to show the number
		 *
		 *  @dtopt Callbacks
		 *  @name DataTable.defaults.formatNumber
		 *
		 *  @example
		 *    // Format a number using a single quote for the separator (note that
		 *    // this can also be done with the language.thousands option)
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "formatNumber": function ( toFormat ) {
		 *          return toFormat.toString().replace(
		 *            /\B(?=(\d{3})+(?!\d))/g, "'"
		 *          );
		 *        };
		 *      } );
		 *    } );
		 */
		"fnFormatNumber": function ( toFormat ) {
			return toFormat.toString().replace(
				/\B(?=(\d{3})+(?!\d))/g,
				this.oLanguage.sThousands
			);
		},
	
	
		/**
		 * This function is called on every 'draw' event, and allows you to
		 * dynamically modify the header row. This can be used to calculate and
		 * display useful information about the table.
		 *  @type function
		 *  @param {node} head "TR" element for the header
		 *  @param {array} data Full table data (as derived from the original HTML)
		 *  @param {int} start Index for the current display starting point in the
		 *    display array
		 *  @param {int} end Index for the current display ending point in the
		 *    display array
		 *  @param {array int} display Index array to translate the visual position
		 *    to the full data array
		 *
		 *  @dtopt Callbacks
		 *  @name DataTable.defaults.headerCallback
		 *
		 *  @example
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "fheaderCallback": function( head, data, start, end, display ) {
		 *          head.getElementsByTagName('th')[0].innerHTML = "Displaying "+(end-start)+" records";
		 *        }
		 *      } );
		 *    } )
		 */
		"fnHeaderCallback": null,
	
	
		/**
		 * The information element can be used to convey information about the current
		 * state of the table. Although the internationalisation options presented by
		 * DataTables are quite capable of dealing with most customisations, there may
		 * be times where you wish to customise the string further. This callback
		 * allows you to do exactly that.
		 *  @type function
		 *  @param {object} oSettings DataTables settings object
		 *  @param {int} start Starting position in data for the draw
		 *  @param {int} end End position in data for the draw
		 *  @param {int} max Total number of rows in the table (regardless of
		 *    filtering)
		 *  @param {int} total Total number of rows in the data set, after filtering
		 *  @param {string} pre The string that DataTables has formatted using it's
		 *    own rules
		 *  @returns {string} The string to be displayed in the information element.
		 *
		 *  @dtopt Callbacks
		 *  @name DataTable.defaults.infoCallback
		 *
		 *  @example
		 *    $('#example').dataTable( {
		 *      "infoCallback": function( settings, start, end, max, total, pre ) {
		 *        return start +" to "+ end;
		 *      }
		 *    } );
		 */
		"fnInfoCallback": null,
	
	
		/**
		 * Called when the table has been initialised. Normally DataTables will
		 * initialise sequentially and there will be no need for this function,
		 * however, this does not hold true when using external language information
		 * since that is obtained using an async XHR call.
		 *  @type function
		 *  @param {object} settings DataTables settings object
		 *  @param {object} json The JSON object request from the server - only
		 *    present if client-side Ajax sourced data is used
		 *
		 *  @dtopt Callbacks
		 *  @name DataTable.defaults.initComplete
		 *
		 *  @example
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "initComplete": function(settings, json) {
		 *          alert( 'DataTables has finished its initialisation.' );
		 *        }
		 *      } );
		 *    } )
		 */
		"fnInitComplete": null,
	
	
		/**
		 * Called at the very start of each table draw and can be used to cancel the
		 * draw by returning false, any other return (including undefined) results in
		 * the full draw occurring).
		 *  @type function
		 *  @param {object} settings DataTables settings object
		 *  @returns {boolean} False will cancel the draw, anything else (including no
		 *    return) will allow it to complete.
		 *
		 *  @dtopt Callbacks
		 *  @name DataTable.defaults.preDrawCallback
		 *
		 *  @example
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "preDrawCallback": function( settings ) {
		 *          if ( $('#test').val() == 1 ) {
		 *            return false;
		 *          }
		 *        }
		 *      } );
		 *    } );
		 */
		"fnPreDrawCallback": null,
	
	
		/**
		 * This function allows you to 'post process' each row after it have been
		 * generated for each table draw, but before it is rendered on screen. This
		 * function might be used for setting the row class name etc.
		 *  @type function
		 *  @param {node} row "TR" element for the current row
		 *  @param {array} data Raw data array for this row
		 *  @param {int} displayIndex The display index for the current table draw
		 *  @param {int} displayIndexFull The index of the data in the full list of
		 *    rows (after filtering)
		 *
		 *  @dtopt Callbacks
		 *  @name DataTable.defaults.rowCallback
		 *
		 *  @example
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "rowCallback": function( row, data, displayIndex, displayIndexFull ) {
		 *          // Bold the grade for all 'A' grade browsers
		 *          if ( data[4] == "A" ) {
		 *            $('td:eq(4)', row).html( '<b>A</b>' );
		 *          }
		 *        }
		 *      } );
		 *    } );
		 */
		"fnRowCallback": null,
	
	
		/**
		 * __Deprecated__ The functionality provided by this parameter has now been
		 * superseded by that provided through `ajax`, which should be used instead.
		 *
		 * This parameter allows you to override the default function which obtains
		 * the data from the server so something more suitable for your application.
		 * For example you could use POST data, or pull information from a Gears or
		 * AIR database.
		 *  @type function
		 *  @member
		 *  @param {string} source HTTP source to obtain the data from (`ajax`)
		 *  @param {array} data A key/value pair object containing the data to send
		 *    to the server
		 *  @param {function} callback to be called on completion of the data get
		 *    process that will draw the data on the page.
		 *  @param {object} settings DataTables settings object
		 *
		 *  @dtopt Callbacks
		 *  @dtopt Server-side
		 *  @name DataTable.defaults.serverData
		 *
		 *  @deprecated 1.10. Please use `ajax` for this functionality now.
		 */
		"fnServerData": null,
	
	
		/**
		 * __Deprecated__ The functionality provided by this parameter has now been
		 * superseded by that provided through `ajax`, which should be used instead.
		 *
		 *  It is often useful to send extra data to the server when making an Ajax
		 * request - for example custom filtering information, and this callback
		 * function makes it trivial to send extra information to the server. The
		 * passed in parameter is the data set that has been constructed by
		 * DataTables, and you can add to this or modify it as you require.
		 *  @type function
		 *  @param {array} data Data array (array of objects which are name/value
		 *    pairs) that has been constructed by DataTables and will be sent to the
		 *    server. In the case of Ajax sourced data with server-side processing
		 *    this will be an empty array, for server-side processing there will be a
		 *    significant number of parameters!
		 *  @returns {undefined} Ensure that you modify the data array passed in,
		 *    as this is passed by reference.
		 *
		 *  @dtopt Callbacks
		 *  @dtopt Server-side
		 *  @name DataTable.defaults.serverParams
		 *
		 *  @deprecated 1.10. Please use `ajax` for this functionality now.
		 */
		"fnServerParams": null,
	
	
		/**
		 * Load the table state. With this function you can define from where, and how, the
		 * state of a table is loaded. By default DataTables will load from `localStorage`
		 * but you might wish to use a server-side database or cookies.
		 *  @type function
		 *  @member
		 *  @param {object} settings DataTables settings object
		 *  @param {object} callback Callback that can be executed when done. It
		 *    should be passed the loaded state object.
		 *  @return {object} The DataTables state object to be loaded
		 *
		 *  @dtopt Callbacks
		 *  @name DataTable.defaults.stateLoadCallback
		 *
		 *  @example
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "stateSave": true,
		 *        "stateLoadCallback": function (settings, callback) {
		 *          $.ajax( {
		 *            "url": "/state_load",
		 *            "dataType": "json",
		 *            "success": function (json) {
		 *              callback( json );
		 *            }
		 *          } );
		 *        }
		 *      } );
		 *    } );
		 */
		"fnStateLoadCallback": function ( settings ) {
			try {
				return JSON.parse(
					(settings.iStateDuration === -1 ? sessionStorage : localStorage).getItem(
						'DataTables_'+settings.sInstance+'_'+location.pathname
					)
				);
			} catch (e) {
				return {};
			}
		},
	
	
		/**
		 * Callback which allows modification of the saved state prior to loading that state.
		 * This callback is called when the table is loading state from the stored data, but
		 * prior to the settings object being modified by the saved state. Note that for
		 * plug-in authors, you should use the `stateLoadParams` event to load parameters for
		 * a plug-in.
		 *  @type function
		 *  @param {object} settings DataTables settings object
		 *  @param {object} data The state object that is to be loaded
		 *
		 *  @dtopt Callbacks
		 *  @name DataTable.defaults.stateLoadParams
		 *
		 *  @example
		 *    // Remove a saved filter, so filtering is never loaded
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "stateSave": true,
		 *        "stateLoadParams": function (settings, data) {
		 *          data.oSearch.sSearch = "";
		 *        }
		 *      } );
		 *    } );
		 *
		 *  @example
		 *    // Disallow state loading by returning false
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "stateSave": true,
		 *        "stateLoadParams": function (settings, data) {
		 *          return false;
		 *        }
		 *      } );
		 *    } );
		 */
		"fnStateLoadParams": null,
	
	
		/**
		 * Callback that is called when the state has been loaded from the state saving method
		 * and the DataTables settings object has been modified as a result of the loaded state.
		 *  @type function
		 *  @param {object} settings DataTables settings object
		 *  @param {object} data The state object that was loaded
		 *
		 *  @dtopt Callbacks
		 *  @name DataTable.defaults.stateLoaded
		 *
		 *  @example
		 *    // Show an alert with the filtering value that was saved
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "stateSave": true,
		 *        "stateLoaded": function (settings, data) {
		 *          alert( 'Saved filter was: '+data.oSearch.sSearch );
		 *        }
		 *      } );
		 *    } );
		 */
		"fnStateLoaded": null,
	
	
		/**
		 * Save the table state. This function allows you to define where and how the state
		 * information for the table is stored By default DataTables will use `localStorage`
		 * but you might wish to use a server-side database or cookies.
		 *  @type function
		 *  @member
		 *  @param {object} settings DataTables settings object
		 *  @param {object} data The state object to be saved
		 *
		 *  @dtopt Callbacks
		 *  @name DataTable.defaults.stateSaveCallback
		 *
		 *  @example
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "stateSave": true,
		 *        "stateSaveCallback": function (settings, data) {
		 *          // Send an Ajax request to the server with the state object
		 *          $.ajax( {
		 *            "url": "/state_save",
		 *            "data": data,
		 *            "dataType": "json",
		 *            "method": "POST"
		 *            "success": function () {}
		 *          } );
		 *        }
		 *      } );
		 *    } );
		 */
		"fnStateSaveCallback": function ( settings, data ) {
			try {
				(settings.iStateDuration === -1 ? sessionStorage : localStorage).setItem(
					'DataTables_'+settings.sInstance+'_'+location.pathname,
					JSON.stringify( data )
				);
			} catch (e) {}
		},
	
	
		/**
		 * Callback which allows modification of the state to be saved. Called when the table
		 * has changed state a new state save is required. This method allows modification of
		 * the state saving object prior to actually doing the save, including addition or
		 * other state properties or modification. Note that for plug-in authors, you should
		 * use the `stateSaveParams` event to save parameters for a plug-in.
		 *  @type function
		 *  @param {object} settings DataTables settings object
		 *  @param {object} data The state object to be saved
		 *
		 *  @dtopt Callbacks
		 *  @name DataTable.defaults.stateSaveParams
		 *
		 *  @example
		 *    // Remove a saved filter, so filtering is never saved
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "stateSave": true,
		 *        "stateSaveParams": function (settings, data) {
		 *          data.oSearch.sSearch = "";
		 *        }
		 *      } );
		 *    } );
		 */
		"fnStateSaveParams": null,
	
	
		/**
		 * Duration for which the saved state information is considered valid. After this period
		 * has elapsed the state will be returned to the default.
		 * Value is given in seconds.
		 *  @type int
		 *  @default 7200 <i>(2 hours)</i>
		 *
		 *  @dtopt Options
		 *  @name DataTable.defaults.stateDuration
		 *
		 *  @example
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "stateDuration": 60*60*24; // 1 day
		 *      } );
		 *    } )
		 */
		"iStateDuration": 7200,
	
	
		/**
		 * When enabled DataTables will not make a request to the server for the first
		 * page draw - rather it will use the data already on the page (no sorting etc
		 * will be applied to it), thus saving on an XHR at load time. `deferLoading`
		 * is used to indicate that deferred loading is required, but it is also used
		 * to tell DataTables how many records there are in the full table (allowing
		 * the information element and pagination to be displayed correctly). In the case
		 * where a filtering is applied to the table on initial load, this can be
		 * indicated by giving the parameter as an array, where the first element is
		 * the number of records available after filtering and the second element is the
		 * number of records without filtering (allowing the table information element
		 * to be shown correctly).
		 *  @type int | array
		 *  @default null
		 *
		 *  @dtopt Options
		 *  @name DataTable.defaults.deferLoading
		 *
		 *  @example
		 *    // 57 records available in the table, no filtering applied
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "serverSide": true,
		 *        "ajax": "scripts/server_processing.php",
		 *        "deferLoading": 57
		 *      } );
		 *    } );
		 *
		 *  @example
		 *    // 57 records after filtering, 100 without filtering (an initial filter applied)
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "serverSide": true,
		 *        "ajax": "scripts/server_processing.php",
		 *        "deferLoading": [ 57, 100 ],
		 *        "search": {
		 *          "search": "my_filter"
		 *        }
		 *      } );
		 *    } );
		 */
		"iDeferLoading": null,
	
	
		/**
		 * Number of rows to display on a single page when using pagination. If
		 * feature enabled (`lengthChange`) then the end user will be able to override
		 * this to a custom setting using a pop-up menu.
		 *  @type int
		 *  @default 10
		 *
		 *  @dtopt Options
		 *  @name DataTable.defaults.pageLength
		 *
		 *  @example
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "pageLength": 50
		 *      } );
		 *    } )
		 */
		"iDisplayLength": 10,
	
	
		/**
		 * Define the starting point for data display when using DataTables with
		 * pagination. Note that this parameter is the number of records, rather than
		 * the page number, so if you have 10 records per page and want to start on
		 * the third page, it should be "20".
		 *  @type int
		 *  @default 0
		 *
		 *  @dtopt Options
		 *  @name DataTable.defaults.displayStart
		 *
		 *  @example
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "displayStart": 20
		 *      } );
		 *    } )
		 */
		"iDisplayStart": 0,
	
	
		/**
		 * By default DataTables allows keyboard navigation of the table (sorting, paging,
		 * and filtering) by adding a `tabindex` attribute to the required elements. This
		 * allows you to tab through the controls and press the enter key to activate them.
		 * The tabindex is default 0, meaning that the tab follows the flow of the document.
		 * You can overrule this using this parameter if you wish. Use a value of -1 to
		 * disable built-in keyboard navigation.
		 *  @type int
		 *  @default 0
		 *
		 *  @dtopt Options
		 *  @name DataTable.defaults.tabIndex
		 *
		 *  @example
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "tabIndex": 1
		 *      } );
		 *    } );
		 */
		"iTabIndex": 0,
	
	
		/**
		 * Classes that DataTables assigns to the various components and features
		 * that it adds to the HTML table. This allows classes to be configured
		 * during initialisation in addition to through the static
		 * {@link DataTable.ext.oStdClasses} object).
		 *  @namespace
		 *  @name DataTable.defaults.classes
		 */
		"oClasses": {},
	
	
		/**
		 * All strings that DataTables uses in the user interface that it creates
		 * are defined in this object, allowing you to modified them individually or
		 * completely replace them all as required.
		 *  @namespace
		 *  @name DataTable.defaults.language
		 */
		"oLanguage": {
			/**
			 * Strings that are used for WAI-ARIA labels and controls only (these are not
			 * actually visible on the page, but will be read by screenreaders, and thus
			 * must be internationalised as well).
			 *  @namespace
			 *  @name DataTable.defaults.language.aria
			 */
			"oAria": {
				/**
				 * ARIA label that is added to the table headers when the column may be
				 * sorted ascending by activing the column (click or return when focused).
				 * Note that the column header is prefixed to this string.
				 *  @type string
				 *  @default : activate to sort column ascending
				 *
				 *  @dtopt Language
				 *  @name DataTable.defaults.language.aria.sortAscending
				 *
				 *  @example
				 *    $(document).ready( function() {
				 *      $('#example').dataTable( {
				 *        "language": {
				 *          "aria": {
				 *            "sortAscending": " - click/return to sort ascending"
				 *          }
				 *        }
				 *      } );
				 *    } );
				 */
				"sSortAscending": ": activate to sort column ascending",
	
				/**
				 * ARIA label that is added to the table headers when the column may be
				 * sorted descending by activing the column (click or return when focused).
				 * Note that the column header is prefixed to this string.
				 *  @type string
				 *  @default : activate to sort column ascending
				 *
				 *  @dtopt Language
				 *  @name DataTable.defaults.language.aria.sortDescending
				 *
				 *  @example
				 *    $(document).ready( function() {
				 *      $('#example').dataTable( {
				 *        "language": {
				 *          "aria": {
				 *            "sortDescending": " - click/return to sort descending"
				 *          }
				 *        }
				 *      } );
				 *    } );
				 */
				"sSortDescending": ": activate to sort column descending"
			},
	
			/**
			 * Pagination string used by DataTables for the built-in pagination
			 * control types.
			 *  @namespace
			 *  @name DataTable.defaults.language.paginate
			 */
			"oPaginate": {
				/**
				 * Text to use when using the 'full_numbers' type of pagination for the
				 * button to take the user to the first page.
				 *  @type string
				 *  @default First
				 *
				 *  @dtopt Language
				 *  @name DataTable.defaults.language.paginate.first
				 *
				 *  @example
				 *    $(document).ready( function() {
				 *      $('#example').dataTable( {
				 *        "language": {
				 *          "paginate": {
				 *            "first": "First page"
				 *          }
				 *        }
				 *      } );
				 *    } );
				 */
				"sFirst": "First",
	
	
				/**
				 * Text to use when using the 'full_numbers' type of pagination for the
				 * button to take the user to the last page.
				 *  @type string
				 *  @default Last
				 *
				 *  @dtopt Language
				 *  @name DataTable.defaults.language.paginate.last
				 *
				 *  @example
				 *    $(document).ready( function() {
				 *      $('#example').dataTable( {
				 *        "language": {
				 *          "paginate": {
				 *            "last": "Last page"
				 *          }
				 *        }
				 *      } );
				 *    } );
				 */
				"sLast": "Last",
	
	
				/**
				 * Text to use for the 'next' pagination button (to take the user to the
				 * next page).
				 *  @type string
				 *  @default Next
				 *
				 *  @dtopt Language
				 *  @name DataTable.defaults.language.paginate.next
				 *
				 *  @example
				 *    $(document).ready( function() {
				 *      $('#example').dataTable( {
				 *        "language": {
				 *          "paginate": {
				 *            "next": "Next page"
				 *          }
				 *        }
				 *      } );
				 *    } );
				 */
				"sNext": "Next",
	
	
				/**
				 * Text to use for the 'previous' pagination button (to take the user to
				 * the previous page).
				 *  @type string
				 *  @default Previous
				 *
				 *  @dtopt Language
				 *  @name DataTable.defaults.language.paginate.previous
				 *
				 *  @example
				 *    $(document).ready( function() {
				 *      $('#example').dataTable( {
				 *        "language": {
				 *          "paginate": {
				 *            "previous": "Previous page"
				 *          }
				 *        }
				 *      } );
				 *    } );
				 */
				"sPrevious": "Previous"
			},
	
			/**
			 * This string is shown in preference to `zeroRecords` when the table is
			 * empty of data (regardless of filtering). Note that this is an optional
			 * parameter - if it is not given, the value of `zeroRecords` will be used
			 * instead (either the default or given value).
			 *  @type string
			 *  @default No data available in table
			 *
			 *  @dtopt Language
			 *  @name DataTable.defaults.language.emptyTable
			 *
			 *  @example
			 *    $(document).ready( function() {
			 *      $('#example').dataTable( {
			 *        "language": {
			 *          "emptyTable": "No data available in table"
			 *        }
			 *      } );
			 *    } );
			 */
			"sEmptyTable": "No data available in table",
	
	
			/**
			 * This string gives information to the end user about the information
			 * that is current on display on the page. The following tokens can be
			 * used in the string and will be dynamically replaced as the table
			 * display updates. This tokens can be placed anywhere in the string, or
			 * removed as needed by the language requires:
			 *
			 * * `\_START\_` - Display index of the first record on the current page
			 * * `\_END\_` - Display index of the last record on the current page
			 * * `\_TOTAL\_` - Number of records in the table after filtering
			 * * `\_MAX\_` - Number of records in the table without filtering
			 * * `\_PAGE\_` - Current page number
			 * * `\_PAGES\_` - Total number of pages of data in the table
			 *
			 *  @type string
			 *  @default Showing _START_ to _END_ of _TOTAL_ entries
			 *
			 *  @dtopt Language
			 *  @name DataTable.defaults.language.info
			 *
			 *  @example
			 *    $(document).ready( function() {
			 *      $('#example').dataTable( {
			 *        "language": {
			 *          "info": "Showing page _PAGE_ of _PAGES_"
			 *        }
			 *      } );
			 *    } );
			 */
			"sInfo": "Showing _START_ to _END_ of _TOTAL_ entries",
	
	
			/**
			 * Display information string for when the table is empty. Typically the
			 * format of this string should match `info`.
			 *  @type string
			 *  @default Showing 0 to 0 of 0 entries
			 *
			 *  @dtopt Language
			 *  @name DataTable.defaults.language.infoEmpty
			 *
			 *  @example
			 *    $(document).ready( function() {
			 *      $('#example').dataTable( {
			 *        "language": {
			 *          "infoEmpty": "No entries to show"
			 *        }
			 *      } );
			 *    } );
			 */
			"sInfoEmpty": "Showing 0 to 0 of 0 entries",
	
	
			/**
			 * When a user filters the information in a table, this string is appended
			 * to the information (`info`) to give an idea of how strong the filtering
			 * is. The variable _MAX_ is dynamically updated.
			 *  @type string
			 *  @default (filtered from _MAX_ total entries)
			 *
			 *  @dtopt Language
			 *  @name DataTable.defaults.language.infoFiltered
			 *
			 *  @example
			 *    $(document).ready( function() {
			 *      $('#example').dataTable( {
			 *        "language": {
			 *          "infoFiltered": " - filtering from _MAX_ records"
			 *        }
			 *      } );
			 *    } );
			 */
			"sInfoFiltered": "(filtered from _MAX_ total entries)",
	
	
			/**
			 * If can be useful to append extra information to the info string at times,
			 * and this variable does exactly that. This information will be appended to
			 * the `info` (`infoEmpty` and `infoFiltered` in whatever combination they are
			 * being used) at all times.
			 *  @type string
			 *  @default <i>Empty string</i>
			 *
			 *  @dtopt Language
			 *  @name DataTable.defaults.language.infoPostFix
			 *
			 *  @example
			 *    $(document).ready( function() {
			 *      $('#example').dataTable( {
			 *        "language": {
			 *          "infoPostFix": "All records shown are derived from real information."
			 *        }
			 *      } );
			 *    } );
			 */
			"sInfoPostFix": "",
	
	
			/**
			 * This decimal place operator is a little different from the other
			 * language options since DataTables doesn't output floating point
			 * numbers, so it won't ever use this for display of a number. Rather,
			 * what this parameter does is modify the sort methods of the table so
			 * that numbers which are in a format which has a character other than
			 * a period (`.`) as a decimal place will be sorted numerically.
			 *
			 * Note that numbers with different decimal places cannot be shown in
			 * the same table and still be sortable, the table must be consistent.
			 * However, multiple different tables on the page can use different
			 * decimal place characters.
			 *  @type string
			 *  @default 
			 *
			 *  @dtopt Language
			 *  @name DataTable.defaults.language.decimal
			 *
			 *  @example
			 *    $(document).ready( function() {
			 *      $('#example').dataTable( {
			 *        "language": {
			 *          "decimal": ","
			 *          "thousands": "."
			 *        }
			 *      } );
			 *    } );
			 */
			"sDecimal": "",
	
	
			/**
			 * DataTables has a build in number formatter (`formatNumber`) which is
			 * used to format large numbers that are used in the table information.
			 * By default a comma is used, but this can be trivially changed to any
			 * character you wish with this parameter.
			 *  @type string
			 *  @default ,
			 *
			 *  @dtopt Language
			 *  @name DataTable.defaults.language.thousands
			 *
			 *  @example
			 *    $(document).ready( function() {
			 *      $('#example').dataTable( {
			 *        "language": {
			 *          "thousands": "'"
			 *        }
			 *      } );
			 *    } );
			 */
			"sThousands": ",",
	
	
			/**
			 * Detail the action that will be taken when the drop down menu for the
			 * pagination length option is changed. The '_MENU_' variable is replaced
			 * with a default select list of 10, 25, 50 and 100, and can be replaced
			 * with a custom select box if required.
			 *  @type string
			 *  @default Show _MENU_ entries
			 *
			 *  @dtopt Language
			 *  @name DataTable.defaults.language.lengthMenu
			 *
			 *  @example
			 *    // Language change only
			 *    $(document).ready( function() {
			 *      $('#example').dataTable( {
			 *        "language": {
			 *          "lengthMenu": "Display _MENU_ records"
			 *        }
			 *      } );
			 *    } );
			 *
			 *  @example
			 *    // Language and options change
			 *    $(document).ready( function() {
			 *      $('#example').dataTable( {
			 *        "language": {
			 *          "lengthMenu": 'Display <select>'+
			 *            '<option value="10">10</option>'+
			 *            '<option value="20">20</option>'+
			 *            '<option value="30">30</option>'+
			 *            '<option value="40">40</option>'+
			 *            '<option value="50">50</option>'+
			 *            '<option value="-1">All</option>'+
			 *            '</select> records'
			 *        }
			 *      } );
			 *    } );
			 */
			"sLengthMenu": "Show _MENU_ entries",
	
	
			/**
			 * When using Ajax sourced data and during the first draw when DataTables is
			 * gathering the data, this message is shown in an empty row in the table to
			 * indicate to the end user the the data is being loaded. Note that this
			 * parameter is not used when loading data by server-side processing, just
			 * Ajax sourced data with client-side processing.
			 *  @type string
			 *  @default Loading...
			 *
			 *  @dtopt Language
			 *  @name DataTable.defaults.language.loadingRecords
			 *
			 *  @example
			 *    $(document).ready( function() {
			 *      $('#example').dataTable( {
			 *        "language": {
			 *          "loadingRecords": "Please wait - loading..."
			 *        }
			 *      } );
			 *    } );
			 */
			"sLoadingRecords": "Loading...",
	
	
			/**
			 * Text which is displayed when the table is processing a user action
			 * (usually a sort command or similar).
			 *  @type string
			 *
			 *  @dtopt Language
			 *  @name DataTable.defaults.language.processing
			 *
			 *  @example
			 *    $(document).ready( function() {
			 *      $('#example').dataTable( {
			 *        "language": {
			 *          "processing": "DataTables is currently busy"
			 *        }
			 *      } );
			 *    } );
			 */
			"sProcessing": "",
	
	
			/**
			 * Details the actions that will be taken when the user types into the
			 * filtering input text box. The variable "_INPUT_", if used in the string,
			 * is replaced with the HTML text box for the filtering input allowing
			 * control over where it appears in the string. If "_INPUT_" is not given
			 * then the input box is appended to the string automatically.
			 *  @type string
			 *  @default Search:
			 *
			 *  @dtopt Language
			 *  @name DataTable.defaults.language.search
			 *
			 *  @example
			 *    // Input text box will be appended at the end automatically
			 *    $(document).ready( function() {
			 *      $('#example').dataTable( {
			 *        "language": {
			 *          "search": "Filter records:"
			 *        }
			 *      } );
			 *    } );
			 *
			 *  @example
			 *    // Specify where the filter should appear
			 *    $(document).ready( function() {
			 *      $('#example').dataTable( {
			 *        "language": {
			 *          "search": "Apply filter _INPUT_ to table"
			 *        }
			 *      } );
			 *    } );
			 */
			"sSearch": "Search:",
	
	
			/**
			 * Assign a `placeholder` attribute to the search `input` element
			 *  @type string
			 *  @default 
			 *
			 *  @dtopt Language
			 *  @name DataTable.defaults.language.searchPlaceholder
			 */
			"sSearchPlaceholder": "",
	
	
			/**
			 * All of the language information can be stored in a file on the
			 * server-side, which DataTables will look up if this parameter is passed.
			 * It must store the URL of the language file, which is in a JSON format,
			 * and the object has the same properties as the oLanguage object in the
			 * initialiser object (i.e. the above parameters). Please refer to one of
			 * the example language files to see how this works in action.
			 *  @type string
			 *  @default <i>Empty string - i.e. disabled</i>
			 *
			 *  @dtopt Language
			 *  @name DataTable.defaults.language.url
			 *
			 *  @example
			 *    $(document).ready( function() {
			 *      $('#example').dataTable( {
			 *        "language": {
			 *          "url": "http://www.sprymedia.co.uk/dataTables/lang.txt"
			 *        }
			 *      } );
			 *    } );
			 */
			"sUrl": "",
	
	
			/**
			 * Text shown inside the table records when the is no information to be
			 * displayed after filtering. `emptyTable` is shown when there is simply no
			 * information in the table at all (regardless of filtering).
			 *  @type string
			 *  @default No matching records found
			 *
			 *  @dtopt Language
			 *  @name DataTable.defaults.language.zeroRecords
			 *
			 *  @example
			 *    $(document).ready( function() {
			 *      $('#example').dataTable( {
			 *        "language": {
			 *          "zeroRecords": "No records to display"
			 *        }
			 *      } );
			 *    } );
			 */
			"sZeroRecords": "No matching records found"
		},
	
	
		/**
		 * This parameter allows you to have define the global filtering state at
		 * initialisation time. As an object the `search` parameter must be
		 * defined, but all other parameters are optional. When `regex` is true,
		 * the search string will be treated as a regular expression, when false
		 * (default) it will be treated as a straight string. When `smart`
		 * DataTables will use it's smart filtering methods (to word match at
		 * any point in the data), when false this will not be done.
		 *  @namespace
		 *  @extends DataTable.models.oSearch
		 *
		 *  @dtopt Options
		 *  @name DataTable.defaults.search
		 *
		 *  @example
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "search": {"search": "Initial search"}
		 *      } );
		 *    } )
		 */
		"oSearch": $.extend( {}, DataTable.models.oSearch ),
	
	
		/**
		 * __Deprecated__ The functionality provided by this parameter has now been
		 * superseded by that provided through `ajax`, which should be used instead.
		 *
		 * By default DataTables will look for the property `data` (or `aaData` for
		 * compatibility with DataTables 1.9-) when obtaining data from an Ajax
		 * source or for server-side processing - this parameter allows that
		 * property to be changed. You can use Javascript dotted object notation to
		 * get a data source for multiple levels of nesting.
		 *  @type string
		 *  @default data
		 *
		 *  @dtopt Options
		 *  @dtopt Server-side
		 *  @name DataTable.defaults.ajaxDataProp
		 *
		 *  @deprecated 1.10. Please use `ajax` for this functionality now.
		 */
		"sAjaxDataProp": "data",
	
	
		/**
		 * __Deprecated__ The functionality provided by this parameter has now been
		 * superseded by that provided through `ajax`, which should be used instead.
		 *
		 * You can instruct DataTables to load data from an external
		 * source using this parameter (use aData if you want to pass data in you
		 * already have). Simply provide a url a JSON object can be obtained from.
		 *  @type string
		 *  @default null
		 *
		 *  @dtopt Options
		 *  @dtopt Server-side
		 *  @name DataTable.defaults.ajaxSource
		 *
		 *  @deprecated 1.10. Please use `ajax` for this functionality now.
		 */
		"sAjaxSource": null,
	
	
		/**
		 * This initialisation variable allows you to specify exactly where in the
		 * DOM you want DataTables to inject the various controls it adds to the page
		 * (for example you might want the pagination controls at the top of the
		 * table). DIV elements (with or without a custom class) can also be added to
		 * aid styling. The follow syntax is used:
		 *   <ul>
		 *     <li>The following options are allowed:
		 *       <ul>
		 *         <li>'l' - Length changing</li>
		 *         <li>'f' - Filtering input</li>
		 *         <li>'t' - The table!</li>
		 *         <li>'i' - Information</li>
		 *         <li>'p' - Pagination</li>
		 *         <li>'r' - pRocessing</li>
		 *       </ul>
		 *     </li>
		 *     <li>The following constants are allowed:
		 *       <ul>
		 *         <li>'H' - jQueryUI theme "header" classes ('fg-toolbar ui-widget-header ui-corner-tl ui-corner-tr ui-helper-clearfix')</li>
		 *         <li>'F' - jQueryUI theme "footer" classes ('fg-toolbar ui-widget-header ui-corner-bl ui-corner-br ui-helper-clearfix')</li>
		 *       </ul>
		 *     </li>
		 *     <li>The following syntax is expected:
		 *       <ul>
		 *         <li>'&lt;' and '&gt;' - div elements</li>
		 *         <li>'&lt;"class" and '&gt;' - div with a class</li>
		 *         <li>'&lt;"#id" and '&gt;' - div with an ID</li>
		 *       </ul>
		 *     </li>
		 *     <li>Examples:
		 *       <ul>
		 *         <li>'&lt;"wrapper"flipt&gt;'</li>
		 *         <li>'&lt;lf&lt;t&gt;ip&gt;'</li>
		 *       </ul>
		 *     </li>
		 *   </ul>
		 *  @type string
		 *  @default lfrtip <i>(when `jQueryUI` is false)</i> <b>or</b>
		 *    <"H"lfr>t<"F"ip> <i>(when `jQueryUI` is true)</i>
		 *
		 *  @dtopt Options
		 *  @name DataTable.defaults.dom
		 *
		 *  @example
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "dom": '&lt;"top"i&gt;rt&lt;"bottom"flp&gt;&lt;"clear"&gt;'
		 *      } );
		 *    } );
		 */
		"sDom": "lfrtip",
	
	
		/**
		 * Search delay option. This will throttle full table searches that use the
		 * DataTables provided search input element (it does not effect calls to
		 * `dt-api search()`, providing a delay before the search is made.
		 *  @type integer
		 *  @default 0
		 *
		 *  @dtopt Options
		 *  @name DataTable.defaults.searchDelay
		 *
		 *  @example
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "searchDelay": 200
		 *      } );
		 *    } )
		 */
		"searchDelay": null,
	
	
		/**
		 * DataTables features six different built-in options for the buttons to
		 * display for pagination control:
		 *
		 * * `numbers` - Page number buttons only
		 * * `simple` - 'Previous' and 'Next' buttons only
		 * * 'simple_numbers` - 'Previous' and 'Next' buttons, plus page numbers
		 * * `full` - 'First', 'Previous', 'Next' and 'Last' buttons
		 * * `full_numbers` - 'First', 'Previous', 'Next' and 'Last' buttons, plus page numbers
		 * * `first_last_numbers` - 'First' and 'Last' buttons, plus page numbers
		 *  
		 * Further methods can be added using {@link DataTable.ext.oPagination}.
		 *  @type string
		 *  @default simple_numbers
		 *
		 *  @dtopt Options
		 *  @name DataTable.defaults.pagingType
		 *
		 *  @example
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "pagingType": "full_numbers"
		 *      } );
		 *    } )
		 */
		"sPaginationType": "simple_numbers",
	
	
		/**
		 * Enable horizontal scrolling. When a table is too wide to fit into a
		 * certain layout, or you have a large number of columns in the table, you
		 * can enable x-scrolling to show the table in a viewport, which can be
		 * scrolled. This property can be `true` which will allow the table to
		 * scroll horizontally when needed, or any CSS unit, or a number (in which
		 * case it will be treated as a pixel measurement). Setting as simply `true`
		 * is recommended.
		 *  @type boolean|string
		 *  @default <i>blank string - i.e. disabled</i>
		 *
		 *  @dtopt Features
		 *  @name DataTable.defaults.scrollX
		 *
		 *  @example
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "scrollX": true,
		 *        "scrollCollapse": true
		 *      } );
		 *    } );
		 */
		"sScrollX": "",
	
	
		/**
		 * This property can be used to force a DataTable to use more width than it
		 * might otherwise do when x-scrolling is enabled. For example if you have a
		 * table which requires to be well spaced, this parameter is useful for
		 * "over-sizing" the table, and thus forcing scrolling. This property can by
		 * any CSS unit, or a number (in which case it will be treated as a pixel
		 * measurement).
		 *  @type string
		 *  @default <i>blank string - i.e. disabled</i>
		 *
		 *  @dtopt Options
		 *  @name DataTable.defaults.scrollXInner
		 *
		 *  @example
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "scrollX": "100%",
		 *        "scrollXInner": "110%"
		 *      } );
		 *    } );
		 */
		"sScrollXInner": "",
	
	
		/**
		 * Enable vertical scrolling. Vertical scrolling will constrain the DataTable
		 * to the given height, and enable scrolling for any data which overflows the
		 * current viewport. This can be used as an alternative to paging to display
		 * a lot of data in a small area (although paging and scrolling can both be
		 * enabled at the same time). This property can be any CSS unit, or a number
		 * (in which case it will be treated as a pixel measurement).
		 *  @type string
		 *  @default <i>blank string - i.e. disabled</i>
		 *
		 *  @dtopt Features
		 *  @name DataTable.defaults.scrollY
		 *
		 *  @example
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "scrollY": "200px",
		 *        "paginate": false
		 *      } );
		 *    } );
		 */
		"sScrollY": "",
	
	
		/**
		 * __Deprecated__ The functionality provided by this parameter has now been
		 * superseded by that provided through `ajax`, which should be used instead.
		 *
		 * Set the HTTP method that is used to make the Ajax call for server-side
		 * processing or Ajax sourced data.
		 *  @type string
		 *  @default GET
		 *
		 *  @dtopt Options
		 *  @dtopt Server-side
		 *  @name DataTable.defaults.serverMethod
		 *
		 *  @deprecated 1.10. Please use `ajax` for this functionality now.
		 */
		"sServerMethod": "GET",
	
	
		/**
		 * DataTables makes use of renderers when displaying HTML elements for
		 * a table. These renderers can be added or modified by plug-ins to
		 * generate suitable mark-up for a site. For example the Bootstrap
		 * integration plug-in for DataTables uses a paging button renderer to
		 * display pagination buttons in the mark-up required by Bootstrap.
		 *
		 * For further information about the renderers available see
		 * DataTable.ext.renderer
		 *  @type string|object
		 *  @default null
		 *
		 *  @name DataTable.defaults.renderer
		 *
		 */
		"renderer": null,
	
	
		/**
		 * Set the data property name that DataTables should use to get a row's id
		 * to set as the `id` property in the node.
		 *  @type string
		 *  @default DT_RowId
		 *
		 *  @name DataTable.defaults.rowId
		 */
		"rowId": "DT_RowId"
	};
	
	_fnHungarianMap( DataTable.defaults );
	
	
	
	/*
	 * Developer note - See note in model.defaults.js about the use of Hungarian
	 * notation and camel case.
	 */
	
	/**
	 * Column options that can be given to DataTables at initialisation time.
	 *  @namespace
	 */
	DataTable.defaults.column = {
		/**
		 * Define which column(s) an order will occur on for this column. This
		 * allows a column's ordering to take multiple columns into account when
		 * doing a sort or use the data from a different column. For example first
		 * name / last name columns make sense to do a multi-column sort over the
		 * two columns.
		 *  @type array|int
		 *  @default null <i>Takes the value of the column index automatically</i>
		 *
		 *  @name DataTable.defaults.column.orderData
		 *  @dtopt Columns
		 *
		 *  @example
		 *    // Using `columnDefs`
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "columnDefs": [
		 *          { "orderData": [ 0, 1 ], "targets": [ 0 ] },
		 *          { "orderData": [ 1, 0 ], "targets": [ 1 ] },
		 *          { "orderData": 2, "targets": [ 2 ] }
		 *        ]
		 *      } );
		 *    } );
		 *
		 *  @example
		 *    // Using `columns`
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "columns": [
		 *          { "orderData": [ 0, 1 ] },
		 *          { "orderData": [ 1, 0 ] },
		 *          { "orderData": 2 },
		 *          null,
		 *          null
		 *        ]
		 *      } );
		 *    } );
		 */
		"aDataSort": null,
		"iDataSort": -1,
	
	
		/**
		 * You can control the default ordering direction, and even alter the
		 * behaviour of the sort handler (i.e. only allow ascending ordering etc)
		 * using this parameter.
		 *  @type array
		 *  @default [ 'asc', 'desc' ]
		 *
		 *  @name DataTable.defaults.column.orderSequence
		 *  @dtopt Columns
		 *
		 *  @example
		 *    // Using `columnDefs`
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "columnDefs": [
		 *          { "orderSequence": [ "asc" ], "targets": [ 1 ] },
		 *          { "orderSequence": [ "desc", "asc", "asc" ], "targets": [ 2 ] },
		 *          { "orderSequence": [ "desc" ], "targets": [ 3 ] }
		 *        ]
		 *      } );
		 *    } );
		 *
		 *  @example
		 *    // Using `columns`
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "columns": [
		 *          null,
		 *          { "orderSequence": [ "asc" ] },
		 *          { "orderSequence": [ "desc", "asc", "asc" ] },
		 *          { "orderSequence": [ "desc" ] },
		 *          null
		 *        ]
		 *      } );
		 *    } );
		 */
		"asSorting": [ 'asc', 'desc' ],
	
	
		/**
		 * Enable or disable filtering on the data in this column.
		 *  @type boolean
		 *  @default true
		 *
		 *  @name DataTable.defaults.column.searchable
		 *  @dtopt Columns
		 *
		 *  @example
		 *    // Using `columnDefs`
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "columnDefs": [
		 *          { "searchable": false, "targets": [ 0 ] }
		 *        ] } );
		 *    } );
		 *
		 *  @example
		 *    // Using `columns`
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "columns": [
		 *          { "searchable": false },
		 *          null,
		 *          null,
		 *          null,
		 *          null
		 *        ] } );
		 *    } );
		 */
		"bSearchable": true,
	
	
		/**
		 * Enable or disable ordering on this column.
		 *  @type boolean
		 *  @default true
		 *
		 *  @name DataTable.defaults.column.orderable
		 *  @dtopt Columns
		 *
		 *  @example
		 *    // Using `columnDefs`
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "columnDefs": [
		 *          { "orderable": false, "targets": [ 0 ] }
		 *        ] } );
		 *    } );
		 *
		 *  @example
		 *    // Using `columns`
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "columns": [
		 *          { "orderable": false },
		 *          null,
		 *          null,
		 *          null,
		 *          null
		 *        ] } );
		 *    } );
		 */
		"bSortable": true,
	
	
		/**
		 * Enable or disable the display of this column.
		 *  @type boolean
		 *  @default true
		 *
		 *  @name DataTable.defaults.column.visible
		 *  @dtopt Columns
		 *
		 *  @example
		 *    // Using `columnDefs`
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "columnDefs": [
		 *          { "visible": false, "targets": [ 0 ] }
		 *        ] } );
		 *    } );
		 *
		 *  @example
		 *    // Using `columns`
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "columns": [
		 *          { "visible": false },
		 *          null,
		 *          null,
		 *          null,
		 *          null
		 *        ] } );
		 *    } );
		 */
		"bVisible": true,
	
	
		/**
		 * Developer definable function that is called whenever a cell is created (Ajax source,
		 * etc) or processed for input (DOM source). This can be used as a compliment to mRender
		 * allowing you to modify the DOM element (add background colour for example) when the
		 * element is available.
		 *  @type function
		 *  @param {element} td The TD node that has been created
		 *  @param {*} cellData The Data for the cell
		 *  @param {array|object} rowData The data for the whole row
		 *  @param {int} row The row index for the aoData data store
		 *  @param {int} col The column index for aoColumns
		 *
		 *  @name DataTable.defaults.column.createdCell
		 *  @dtopt Columns
		 *
		 *  @example
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "columnDefs": [ {
		 *          "targets": [3],
		 *          "createdCell": function (td, cellData, rowData, row, col) {
		 *            if ( cellData == "1.7" ) {
		 *              $(td).css('color', 'blue')
		 *            }
		 *          }
		 *        } ]
		 *      });
		 *    } );
		 */
		"fnCreatedCell": null,
	
	
		/**
		 * This parameter has been replaced by `data` in DataTables to ensure naming
		 * consistency. `dataProp` can still be used, as there is backwards
		 * compatibility in DataTables for this option, but it is strongly
		 * recommended that you use `data` in preference to `dataProp`.
		 *  @name DataTable.defaults.column.dataProp
		 */
	
	
		/**
		 * This property can be used to read data from any data source property,
		 * including deeply nested objects / properties. `data` can be given in a
		 * number of different ways which effect its behaviour:
		 *
		 * * `integer` - treated as an array index for the data source. This is the
		 *   default that DataTables uses (incrementally increased for each column).
		 * * `string` - read an object property from the data source. There are
		 *   three 'special' options that can be used in the string to alter how
		 *   DataTables reads the data from the source object:
		 *    * `.` - Dotted Javascript notation. Just as you use a `.` in
		 *      Javascript to read from nested objects, so to can the options
		 *      specified in `data`. For example: `browser.version` or
		 *      `browser.name`. If your object parameter name contains a period, use
		 *      `\\` to escape it - i.e. `first\\.name`.
		 *    * `[]` - Array notation. DataTables can automatically combine data
		 *      from and array source, joining the data with the characters provided
		 *      between the two brackets. For example: `name[, ]` would provide a
		 *      comma-space separated list from the source array. If no characters
		 *      are provided between the brackets, the original array source is
		 *      returned.
		 *    * `()` - Function notation. Adding `()` to the end of a parameter will
		 *      execute a function of the name given. For example: `browser()` for a
		 *      simple function on the data source, `browser.version()` for a
		 *      function in a nested property or even `browser().version` to get an
		 *      object property if the function called returns an object. Note that
		 *      function notation is recommended for use in `render` rather than
		 *      `data` as it is much simpler to use as a renderer.
		 * * `null` - use the original data source for the row rather than plucking
		 *   data directly from it. This action has effects on two other
		 *   initialisation options:
		 *    * `defaultContent` - When null is given as the `data` option and
		 *      `defaultContent` is specified for the column, the value defined by
		 *      `defaultContent` will be used for the cell.
		 *    * `render` - When null is used for the `data` option and the `render`
		 *      option is specified for the column, the whole data source for the
		 *      row is used for the renderer.
		 * * `function` - the function given will be executed whenever DataTables
		 *   needs to set or get the data for a cell in the column. The function
		 *   takes three parameters:
		 *    * Parameters:
		 *      * `{array|object}` The data source for the row
		 *      * `{string}` The type call data requested - this will be 'set' when
		 *        setting data or 'filter', 'display', 'type', 'sort' or undefined
		 *        when gathering data. Note that when `undefined` is given for the
		 *        type DataTables expects to get the raw data for the object back<
		 *      * `{*}` Data to set when the second parameter is 'set'.
		 *    * Return:
		 *      * The return value from the function is not required when 'set' is
		 *        the type of call, but otherwise the return is what will be used
		 *        for the data requested.
		 *
		 * Note that `data` is a getter and setter option. If you just require
		 * formatting of data for output, you will likely want to use `render` which
		 * is simply a getter and thus simpler to use.
		 *
		 * Note that prior to DataTables 1.9.2 `data` was called `mDataProp`. The
		 * name change reflects the flexibility of this property and is consistent
		 * with the naming of mRender. If 'mDataProp' is given, then it will still
		 * be used by DataTables, as it automatically maps the old name to the new
		 * if required.
		 *
		 *  @type string|int|function|null
		 *  @default null <i>Use automatically calculated column index</i>
		 *
		 *  @name DataTable.defaults.column.data
		 *  @dtopt Columns
		 *
		 *  @example
		 *    // Read table data from objects
		 *    // JSON structure for each row:
		 *    //   {
		 *    //      "engine": {value},
		 *    //      "browser": {value},
		 *    //      "platform": {value},
		 *    //      "version": {value},
		 *    //      "grade": {value}
		 *    //   }
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "ajaxSource": "sources/objects.txt",
		 *        "columns": [
		 *          { "data": "engine" },
		 *          { "data": "browser" },
		 *          { "data": "platform" },
		 *          { "data": "version" },
		 *          { "data": "grade" }
		 *        ]
		 *      } );
		 *    } );
		 *
		 *  @example
		 *    // Read information from deeply nested objects
		 *    // JSON structure for each row:
		 *    //   {
		 *    //      "engine": {value},
		 *    //      "browser": {value},
		 *    //      "platform": {
		 *    //         "inner": {value}
		 *    //      },
		 *    //      "details": [
		 *    //         {value}, {value}
		 *    //      ]
		 *    //   }
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "ajaxSource": "sources/deep.txt",
		 *        "columns": [
		 *          { "data": "engine" },
		 *          { "data": "browser" },
		 *          { "data": "platform.inner" },
		 *          { "data": "details.0" },
		 *          { "data": "details.1" }
		 *        ]
		 *      } );
		 *    } );
		 *
		 *  @example
		 *    // Using `data` as a function to provide different information for
		 *    // sorting, filtering and display. In this case, currency (price)
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "columnDefs": [ {
		 *          "targets": [ 0 ],
		 *          "data": function ( source, type, val ) {
		 *            if (type === 'set') {
		 *              source.price = val;
		 *              // Store the computed display and filter values for efficiency
		 *              source.price_display = val=="" ? "" : "$"+numberFormat(val);
		 *              source.price_filter  = val=="" ? "" : "$"+numberFormat(val)+" "+val;
		 *              return;
		 *            }
		 *            else if (type === 'display') {
		 *              return source.price_display;
		 *            }
		 *            else if (type === 'filter') {
		 *              return source.price_filter;
		 *            }
		 *            // 'sort', 'type' and undefined all just use the integer
		 *            return source.price;
		 *          }
		 *        } ]
		 *      } );
		 *    } );
		 *
		 *  @example
		 *    // Using default content
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "columnDefs": [ {
		 *          "targets": [ 0 ],
		 *          "data": null,
		 *          "defaultContent": "Click to edit"
		 *        } ]
		 *      } );
		 *    } );
		 *
		 *  @example
		 *    // Using array notation - outputting a list from an array
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "columnDefs": [ {
		 *          "targets": [ 0 ],
		 *          "data": "name[, ]"
		 *        } ]
		 *      } );
		 *    } );
		 *
		 */
		"mData": null,
	
	
		/**
		 * This property is the rendering partner to `data` and it is suggested that
		 * when you want to manipulate data for display (including filtering,
		 * sorting etc) without altering the underlying data for the table, use this
		 * property. `render` can be considered to be the the read only companion to
		 * `data` which is read / write (then as such more complex). Like `data`
		 * this option can be given in a number of different ways to effect its
		 * behaviour:
		 *
		 * * `integer` - treated as an array index for the data source. This is the
		 *   default that DataTables uses (incrementally increased for each column).
		 * * `string` - read an object property from the data source. There are
		 *   three 'special' options that can be used in the string to alter how
		 *   DataTables reads the data from the source object:
		 *    * `.` - Dotted Javascript notation. Just as you use a `.` in
		 *      Javascript to read from nested objects, so to can the options
		 *      specified in `data`. For example: `browser.version` or
		 *      `browser.name`. If your object parameter name contains a period, use
		 *      `\\` to escape it - i.e. `first\\.name`.
		 *    * `[]` - Array notation. DataTables can automatically combine data
		 *      from and array source, joining the data with the characters provided
		 *      between the two brackets. For example: `name[, ]` would provide a
		 *      comma-space separated list from the source array. If no characters
		 *      are provided between the brackets, the original array source is
		 *      returned.
		 *    * `()` - Function notation. Adding `()` to the end of a parameter will
		 *      execute a function of the name given. For example: `browser()` for a
		 *      simple function on the data source, `browser.version()` for a
		 *      function in a nested property or even `browser().version` to get an
		 *      object property if the function called returns an object.
		 * * `object` - use different data for the different data types requested by
		 *   DataTables ('filter', 'display', 'type' or 'sort'). The property names
		 *   of the object is the data type the property refers to and the value can
		 *   defined using an integer, string or function using the same rules as
		 *   `render` normally does. Note that an `_` option _must_ be specified.
		 *   This is the default value to use if you haven't specified a value for
		 *   the data type requested by DataTables.
		 * * `function` - the function given will be executed whenever DataTables
		 *   needs to set or get the data for a cell in the column. The function
		 *   takes three parameters:
		 *    * Parameters:
		 *      * {array|object} The data source for the row (based on `data`)
		 *      * {string} The type call data requested - this will be 'filter',
		 *        'display', 'type' or 'sort'.
		 *      * {array|object} The full data source for the row (not based on
		 *        `data`)
		 *    * Return:
		 *      * The return value from the function is what will be used for the
		 *        data requested.
		 *
		 *  @type string|int|function|object|null
		 *  @default null Use the data source value.
		 *
		 *  @name DataTable.defaults.column.render
		 *  @dtopt Columns
		 *
		 *  @example
		 *    // Create a comma separated list from an array of objects
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "ajaxSource": "sources/deep.txt",
		 *        "columns": [
		 *          { "data": "engine" },
		 *          { "data": "browser" },
		 *          {
		 *            "data": "platform",
		 *            "render": "[, ].name"
		 *          }
		 *        ]
		 *      } );
		 *    } );
		 *
		 *  @example
		 *    // Execute a function to obtain data
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "columnDefs": [ {
		 *          "targets": [ 0 ],
		 *          "data": null, // Use the full data source object for the renderer's source
		 *          "render": "browserName()"
		 *        } ]
		 *      } );
		 *    } );
		 *
		 *  @example
		 *    // As an object, extracting different data for the different types
		 *    // This would be used with a data source such as:
		 *    //   { "phone": 5552368, "phone_filter": "5552368 555-2368", "phone_display": "555-2368" }
		 *    // Here the `phone` integer is used for sorting and type detection, while `phone_filter`
		 *    // (which has both forms) is used for filtering for if a user inputs either format, while
		 *    // the formatted phone number is the one that is shown in the table.
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "columnDefs": [ {
		 *          "targets": [ 0 ],
		 *          "data": null, // Use the full data source object for the renderer's source
		 *          "render": {
		 *            "_": "phone",
		 *            "filter": "phone_filter",
		 *            "display": "phone_display"
		 *          }
		 *        } ]
		 *      } );
		 *    } );
		 *
		 *  @example
		 *    // Use as a function to create a link from the data source
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "columnDefs": [ {
		 *          "targets": [ 0 ],
		 *          "data": "download_link",
		 *          "render": function ( data, type, full ) {
		 *            return '<a href="'+data+'">Download</a>';
		 *          }
		 *        } ]
		 *      } );
		 *    } );
		 */
		"mRender": null,
	
	
		/**
		 * Change the cell type created for the column - either TD cells or TH cells. This
		 * can be useful as TH cells have semantic meaning in the table body, allowing them
		 * to act as a header for a row (you may wish to add scope='row' to the TH elements).
		 *  @type string
		 *  @default td
		 *
		 *  @name DataTable.defaults.column.cellType
		 *  @dtopt Columns
		 *
		 *  @example
		 *    // Make the first column use TH cells
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "columnDefs": [ {
		 *          "targets": [ 0 ],
		 *          "cellType": "th"
		 *        } ]
		 *      } );
		 *    } );
		 */
		"sCellType": "td",
	
	
		/**
		 * Class to give to each cell in this column.
		 *  @type string
		 *  @default <i>Empty string</i>
		 *
		 *  @name DataTable.defaults.column.class
		 *  @dtopt Columns
		 *
		 *  @example
		 *    // Using `columnDefs`
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "columnDefs": [
		 *          { "class": "my_class", "targets": [ 0 ] }
		 *        ]
		 *      } );
		 *    } );
		 *
		 *  @example
		 *    // Using `columns`
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "columns": [
		 *          { "class": "my_class" },
		 *          null,
		 *          null,
		 *          null,
		 *          null
		 *        ]
		 *      } );
		 *    } );
		 */
		"sClass": "",
	
		/**
		 * When DataTables calculates the column widths to assign to each column,
		 * it finds the longest string in each column and then constructs a
		 * temporary table and reads the widths from that. The problem with this
		 * is that "mmm" is much wider then "iiii", but the latter is a longer
		 * string - thus the calculation can go wrong (doing it properly and putting
		 * it into an DOM object and measuring that is horribly(!) slow). Thus as
		 * a "work around" we provide this option. It will append its value to the
		 * text that is found to be the longest string for the column - i.e. padding.
		 * Generally you shouldn't need this!
		 *  @type string
		 *  @default <i>Empty string<i>
		 *
		 *  @name DataTable.defaults.column.contentPadding
		 *  @dtopt Columns
		 *
		 *  @example
		 *    // Using `columns`
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "columns": [
		 *          null,
		 *          null,
		 *          null,
		 *          {
		 *            "contentPadding": "mmm"
		 *          }
		 *        ]
		 *      } );
		 *    } );
		 */
		"sContentPadding": "",
	
	
		/**
		 * Allows a default value to be given for a column's data, and will be used
		 * whenever a null data source is encountered (this can be because `data`
		 * is set to null, or because the data source itself is null).
		 *  @type string
		 *  @default null
		 *
		 *  @name DataTable.defaults.column.defaultContent
		 *  @dtopt Columns
		 *
		 *  @example
		 *    // Using `columnDefs`
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "columnDefs": [
		 *          {
		 *            "data": null,
		 *            "defaultContent": "Edit",
		 *            "targets": [ -1 ]
		 *          }
		 *        ]
		 *      } );
		 *    } );
		 *
		 *  @example
		 *    // Using `columns`
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "columns": [
		 *          null,
		 *          null,
		 *          null,
		 *          {
		 *            "data": null,
		 *            "defaultContent": "Edit"
		 *          }
		 *        ]
		 *      } );
		 *    } );
		 */
		"sDefaultContent": null,
	
	
		/**
		 * This parameter is only used in DataTables' server-side processing. It can
		 * be exceptionally useful to know what columns are being displayed on the
		 * client side, and to map these to database fields. When defined, the names
		 * also allow DataTables to reorder information from the server if it comes
		 * back in an unexpected order (i.e. if you switch your columns around on the
		 * client-side, your server-side code does not also need updating).
		 *  @type string
		 *  @default <i>Empty string</i>
		 *
		 *  @name DataTable.defaults.column.name
		 *  @dtopt Columns
		 *
		 *  @example
		 *    // Using `columnDefs`
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "columnDefs": [
		 *          { "name": "engine", "targets": [ 0 ] },
		 *          { "name": "browser", "targets": [ 1 ] },
		 *          { "name": "platform", "targets": [ 2 ] },
		 *          { "name": "version", "targets": [ 3 ] },
		 *          { "name": "grade", "targets": [ 4 ] }
		 *        ]
		 *      } );
		 *    } );
		 *
		 *  @example
		 *    // Using `columns`
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "columns": [
		 *          { "name": "engine" },
		 *          { "name": "browser" },
		 *          { "name": "platform" },
		 *          { "name": "version" },
		 *          { "name": "grade" }
		 *        ]
		 *      } );
		 *    } );
		 */
		"sName": "",
	
	
		/**
		 * Defines a data source type for the ordering which can be used to read
		 * real-time information from the table (updating the internally cached
		 * version) prior to ordering. This allows ordering to occur on user
		 * editable elements such as form inputs.
		 *  @type string
		 *  @default std
		 *
		 *  @name DataTable.defaults.column.orderDataType
		 *  @dtopt Columns
		 *
		 *  @example
		 *    // Using `columnDefs`
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "columnDefs": [
		 *          { "orderDataType": "dom-text", "targets": [ 2, 3 ] },
		 *          { "type": "numeric", "targets": [ 3 ] },
		 *          { "orderDataType": "dom-select", "targets": [ 4 ] },
		 *          { "orderDataType": "dom-checkbox", "targets": [ 5 ] }
		 *        ]
		 *      } );
		 *    } );
		 *
		 *  @example
		 *    // Using `columns`
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "columns": [
		 *          null,
		 *          null,
		 *          { "orderDataType": "dom-text" },
		 *          { "orderDataType": "dom-text", "type": "numeric" },
		 *          { "orderDataType": "dom-select" },
		 *          { "orderDataType": "dom-checkbox" }
		 *        ]
		 *      } );
		 *    } );
		 */
		"sSortDataType": "std",
	
	
		/**
		 * The title of this column.
		 *  @type string
		 *  @default null <i>Derived from the 'TH' value for this column in the
		 *    original HTML table.</i>
		 *
		 *  @name DataTable.defaults.column.title
		 *  @dtopt Columns
		 *
		 *  @example
		 *    // Using `columnDefs`
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "columnDefs": [
		 *          { "title": "My column title", "targets": [ 0 ] }
		 *        ]
		 *      } );
		 *    } );
		 *
		 *  @example
		 *    // Using `columns`
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "columns": [
		 *          { "title": "My column title" },
		 *          null,
		 *          null,
		 *          null,
		 *          null
		 *        ]
		 *      } );
		 *    } );
		 */
		"sTitle": null,
	
	
		/**
		 * The type allows you to specify how the data for this column will be
		 * ordered. Four types (string, numeric, date and html (which will strip
		 * HTML tags before ordering)) are currently available. Note that only date
		 * formats understood by Javascript's Date() object will be accepted as type
		 * date. For example: "Mar 26, 2008 5:03 PM". May take the values: 'string',
		 * 'numeric', 'date' or 'html' (by default). Further types can be adding
		 * through plug-ins.
		 *  @type string
		 *  @default null <i>Auto-detected from raw data</i>
		 *
		 *  @name DataTable.defaults.column.type
		 *  @dtopt Columns
		 *
		 *  @example
		 *    // Using `columnDefs`
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "columnDefs": [
		 *          { "type": "html", "targets": [ 0 ] }
		 *        ]
		 *      } );
		 *    } );
		 *
		 *  @example
		 *    // Using `columns`
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "columns": [
		 *          { "type": "html" },
		 *          null,
		 *          null,
		 *          null,
		 *          null
		 *        ]
		 *      } );
		 *    } );
		 */
		"sType": null,
	
	
		/**
		 * Defining the width of the column, this parameter may take any CSS value
		 * (3em, 20px etc). DataTables applies 'smart' widths to columns which have not
		 * been given a specific width through this interface ensuring that the table
		 * remains readable.
		 *  @type string
		 *  @default null <i>Automatic</i>
		 *
		 *  @name DataTable.defaults.column.width
		 *  @dtopt Columns
		 *
		 *  @example
		 *    // Using `columnDefs`
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "columnDefs": [
		 *          { "width": "20%", "targets": [ 0 ] }
		 *        ]
		 *      } );
		 *    } );
		 *
		 *  @example
		 *    // Using `columns`
		 *    $(document).ready( function() {
		 *      $('#example').dataTable( {
		 *        "columns": [
		 *          { "width": "20%" },
		 *          null,
		 *          null,
		 *          null,
		 *          null
		 *        ]
		 *      } );
		 *    } );
		 */
		"sWidth": null
	};
	
	_fnHungarianMap( DataTable.defaults.column );
	
	
	
	/**
	 * DataTables settings object - this holds all the information needed for a
	 * given table, including configuration, data and current application of the
	 * table options. DataTables does not have a single instance for each DataTable
	 * with the settings attached to that instance, but rather instances of the
	 * DataTable "class" are created on-the-fly as needed (typically by a
	 * $().dataTable() call) and the settings object is then applied to that
	 * instance.
	 *
	 * Note that this object is related to {@link DataTable.defaults} but this
	 * one is the internal data store for DataTables's cache of columns. It should
	 * NOT be manipulated outside of DataTables. Any configuration should be done
	 * through the initialisation options.
	 *  @namespace
	 *  @todo Really should attach the settings object to individual instances so we
	 *    don't need to create new instances on each $().dataTable() call (if the
	 *    table already exists). It would also save passing oSettings around and
	 *    into every single function. However, this is a very significant
	 *    architecture change for DataTables and will almost certainly break
	 *    backwards compatibility with older installations. This is something that
	 *    will be done in 2.0.
	 */
	DataTable.models.oSettings = {
		/**
		 * Primary features of DataTables and their enablement state.
		 *  @namespace
		 */
		"oFeatures": {
	
			/**
			 * Flag to say if DataTables should automatically try to calculate the
			 * optimum table and columns widths (true) or not (false).
			 * Note that this parameter will be set by the initialisation routine. To
			 * set a default use {@link DataTable.defaults}.
			 *  @type boolean
			 */
			"bAutoWidth": null,
	
			/**
			 * Delay the creation of TR and TD elements until they are actually
			 * needed by a driven page draw. This can give a significant speed
			 * increase for Ajax source and Javascript source data, but makes no
			 * difference at all for DOM and server-side processing tables.
			 * Note that this parameter will be set by the initialisation routine. To
			 * set a default use {@link DataTable.defaults}.
			 *  @type boolean
			 */
			"bDeferRender": null,
	
			/**
			 * Enable filtering on the table or not. Note that if this is disabled
			 * then there is no filtering at all on the table, including fnFilter.
			 * To just remove the filtering input use sDom and remove the 'f' option.
			 * Note that this parameter will be set by the initialisation routine. To
			 * set a default use {@link DataTable.defaults}.
			 *  @type boolean
			 */
			"bFilter": null,
	
			/**
			 * Table information element (the 'Showing x of y records' div) enable
			 * flag.
			 * Note that this parameter will be set by the initialisation routine. To
			 * set a default use {@link DataTable.defaults}.
			 *  @type boolean
			 */
			"bInfo": null,
	
			/**
			 * Present a user control allowing the end user to change the page size
			 * when pagination is enabled.
			 * Note that this parameter will be set by the initialisation routine. To
			 * set a default use {@link DataTable.defaults}.
			 *  @type boolean
			 */
			"bLengthChange": null,
	
			/**
			 * Pagination enabled or not. Note that if this is disabled then length
			 * changing must also be disabled.
			 * Note that this parameter will be set by the initialisation routine. To
			 * set a default use {@link DataTable.defaults}.
			 *  @type boolean
			 */
			"bPaginate": null,
	
			/**
			 * Processing indicator enable flag whenever DataTables is enacting a
			 * user request - typically an Ajax request for server-side processing.
			 * Note that this parameter will be set by the initialisation routine. To
			 * set a default use {@link DataTable.defaults}.
			 *  @type boolean
			 */
			"bProcessing": null,
	
			/**
			 * Server-side processing enabled flag - when enabled DataTables will
			 * get all data from the server for every draw - there is no filtering,
			 * sorting or paging done on the client-side.
			 * Note that this parameter will be set by the initialisation routine. To
			 * set a default use {@link DataTable.defaults}.
			 *  @type boolean
			 */
			"bServerSide": null,
	
			/**
			 * Sorting enablement flag.
			 * Note that this parameter will be set by the initialisation routine. To
			 * set a default use {@link DataTable.defaults}.
			 *  @type boolean
			 */
			"bSort": null,
	
			/**
			 * Multi-column sorting
			 * Note that this parameter will be set by the initialisation routine. To
			 * set a default use {@link DataTable.defaults}.
			 *  @type boolean
			 */
			"bSortMulti": null,
	
			/**
			 * Apply a class to the columns which are being sorted to provide a
			 * visual highlight or not. This can slow things down when enabled since
			 * there is a lot of DOM interaction.
			 * Note that this parameter will be set by the initialisation routine. To
			 * set a default use {@link DataTable.defaults}.
			 *  @type boolean
			 */
			"bSortClasses": null,
	
			/**
			 * State saving enablement flag.
			 * Note that this parameter will be set by the initialisation routine. To
			 * set a default use {@link DataTable.defaults}.
			 *  @type boolean
			 */
			"bStateSave": null
		},
	
	
		/**
		 * Scrolling settings for a table.
		 *  @namespace
		 */
		"oScroll": {
			/**
			 * When the table is shorter in height than sScrollY, collapse the
			 * table container down to the height of the table (when true).
			 * Note that this parameter will be set by the initialisation routine. To
			 * set a default use {@link DataTable.defaults}.
			 *  @type boolean
			 */
			"bCollapse": null,
	
			/**
			 * Width of the scrollbar for the web-browser's platform. Calculated
			 * during table initialisation.
			 *  @type int
			 *  @default 0
			 */
			"iBarWidth": 0,
	
			/**
			 * Viewport width for horizontal scrolling. Horizontal scrolling is
			 * disabled if an empty string.
			 * Note that this parameter will be set by the initialisation routine. To
			 * set a default use {@link DataTable.defaults}.
			 *  @type string
			 */
			"sX": null,
	
			/**
			 * Width to expand the table to when using x-scrolling. Typically you
			 * should not need to use this.
			 * Note that this parameter will be set by the initialisation routine. To
			 * set a default use {@link DataTable.defaults}.
			 *  @type string
			 *  @deprecated
			 */
			"sXInner": null,
	
			/**
			 * Viewport height for vertical scrolling. Vertical scrolling is disabled
			 * if an empty string.
			 * Note that this parameter will be set by the initialisation routine. To
			 * set a default use {@link DataTable.defaults}.
			 *  @type string
			 */
			"sY": null
		},
	
		/**
		 * Language information for the table.
		 *  @namespace
		 *  @extends DataTable.defaults.oLanguage
		 */
		"oLanguage": {
			/**
			 * Information callback function. See
			 * {@link DataTable.defaults.fnInfoCallback}
			 *  @type function
			 *  @default null
			 */
			"fnInfoCallback": null
		},
	
		/**
		 * Browser support parameters
		 *  @namespace
		 */
		"oBrowser": {
			/**
			 * Indicate if the browser incorrectly calculates width:100% inside a
			 * scrolling element (IE6/7)
			 *  @type boolean
			 *  @default false
			 */
			"bScrollOversize": false,
	
			/**
			 * Determine if the vertical scrollbar is on the right or left of the
			 * scrolling container - needed for rtl language layout, although not
			 * all browsers move the scrollbar (Safari).
			 *  @type boolean
			 *  @default false
			 */
			"bScrollbarLeft": false,
	
			/**
			 * Flag for if `getBoundingClientRect` is fully supported or not
			 *  @type boolean
			 *  @default false
			 */
			"bBounding": false,
	
			/**
			 * Browser scrollbar width
			 *  @type integer
			 *  @default 0
			 */
			"barWidth": 0
		},
	
	
		"ajax": null,
	
	
		/**
		 * Array referencing the nodes which are used for the features. The
		 * parameters of this object match what is allowed by sDom - i.e.
		 *   <ul>
		 *     <li>'l' - Length changing</li>
		 *     <li>'f' - Filtering input</li>
		 *     <li>'t' - The table!</li>
		 *     <li>'i' - Information</li>
		 *     <li>'p' - Pagination</li>
		 *     <li>'r' - pRocessing</li>
		 *   </ul>
		 *  @type array
		 *  @default []
		 */
		"aanFeatures": [],
	
		/**
		 * Store data information - see {@link DataTable.models.oRow} for detailed
		 * information.
		 *  @type array
		 *  @default []
		 */
		"aoData": [],
	
		/**
		 * Array of indexes which are in the current display (after filtering etc)
		 *  @type array
		 *  @default []
		 */
		"aiDisplay": [],
	
		/**
		 * Array of indexes for display - no filtering
		 *  @type array
		 *  @default []
		 */
		"aiDisplayMaster": [],
	
		/**
		 * Map of row ids to data indexes
		 *  @type object
		 *  @default {}
		 */
		"aIds": {},
	
		/**
		 * Store information about each column that is in use
		 *  @type array
		 *  @default []
		 */
		"aoColumns": [],
	
		/**
		 * Store information about the table's header
		 *  @type array
		 *  @default []
		 */
		"aoHeader": [],
	
		/**
		 * Store information about the table's footer
		 *  @type array
		 *  @default []
		 */
		"aoFooter": [],
	
		/**
		 * Store the applied global search information in case we want to force a
		 * research or compare the old search to a new one.
		 * Note that this parameter will be set by the initialisation routine. To
		 * set a default use {@link DataTable.defaults}.
		 *  @namespace
		 *  @extends DataTable.models.oSearch
		 */
		"oPreviousSearch": {},
	
		/**
		 * Store the applied search for each column - see
		 * {@link DataTable.models.oSearch} for the format that is used for the
		 * filtering information for each column.
		 *  @type array
		 *  @default []
		 */
		"aoPreSearchCols": [],
	
		/**
		 * Sorting that is applied to the table. Note that the inner arrays are
		 * used in the following manner:
		 * <ul>
		 *   <li>Index 0 - column number</li>
		 *   <li>Index 1 - current sorting direction</li>
		 * </ul>
		 * Note that this parameter will be set by the initialisation routine. To
		 * set a default use {@link DataTable.defaults}.
		 *  @type array
		 *  @todo These inner arrays should really be objects
		 */
		"aaSorting": null,
	
		/**
		 * Sorting that is always applied to the table (i.e. prefixed in front of
		 * aaSorting).
		 * Note that this parameter will be set by the initialisation routine. To
		 * set a default use {@link DataTable.defaults}.
		 *  @type array
		 *  @default []
		 */
		"aaSortingFixed": [],
	
		/**
		 * Classes to use for the striping of a table.
		 * Note that this parameter will be set by the initialisation routine. To
		 * set a default use {@link DataTable.defaults}.
		 *  @type array
		 *  @default []
		 */
		"asStripeClasses": null,
	
		/**
		 * If restoring a table - we should restore its striping classes as well
		 *  @type array
		 *  @default []
		 */
		"asDestroyStripes": [],
	
		/**
		 * If restoring a table - we should restore its width
		 *  @type int
		 *  @default 0
		 */
		"sDestroyWidth": 0,
	
		/**
		 * Callback functions array for every time a row is inserted (i.e. on a draw).
		 *  @type array
		 *  @default []
		 */
		"aoRowCallback": [],
	
		/**
		 * Callback functions for the header on each draw.
		 *  @type array
		 *  @default []
		 */
		"aoHeaderCallback": [],
	
		/**
		 * Callback function for the footer on each draw.
		 *  @type array
		 *  @default []
		 */
		"aoFooterCallback": [],
	
		/**
		 * Array of callback functions for draw callback functions
		 *  @type array
		 *  @default []
		 */
		"aoDrawCallback": [],
	
		/**
		 * Array of callback functions for row created function
		 *  @type array
		 *  @default []
		 */
		"aoRowCreatedCallback": [],
	
		/**
		 * Callback functions for just before the table is redrawn. A return of
		 * false will be used to cancel the draw.
		 *  @type array
		 *  @default []
		 */
		"aoPreDrawCallback": [],
	
		/**
		 * Callback functions for when the table has been initialised.
		 *  @type array
		 *  @default []
		 */
		"aoInitComplete": [],
	
	
		/**
		 * Callbacks for modifying the settings to be stored for state saving, prior to
		 * saving state.
		 *  @type array
		 *  @default []
		 */
		"aoStateSaveParams": [],
	
		/**
		 * Callbacks for modifying the settings that have been stored for state saving
		 * prior to using the stored values to restore the state.
		 *  @type array
		 *  @default []
		 */
		"aoStateLoadParams": [],
	
		/**
		 * Callbacks for operating on the settings object once the saved state has been
		 * loaded
		 *  @type array
		 *  @default []
		 */
		"aoStateLoaded": [],
	
		/**
		 * Cache the table ID for quick access
		 *  @type string
		 *  @default <i>Empty string</i>
		 */
		"sTableId": "",
	
		/**
		 * The TABLE node for the main table
		 *  @type node
		 *  @default null
		 */
		"nTable": null,
	
		/**
		 * Permanent ref to the thead element
		 *  @type node
		 *  @default null
		 */
		"nTHead": null,
	
		/**
		 * Permanent ref to the tfoot element - if it exists
		 *  @type node
		 *  @default null
		 */
		"nTFoot": null,
	
		/**
		 * Permanent ref to the tbody element
		 *  @type node
		 *  @default null
		 */
		"nTBody": null,
	
		/**
		 * Cache the wrapper node (contains all DataTables controlled elements)
		 *  @type node
		 *  @default null
		 */
		"nTableWrapper": null,
	
		/**
		 * Indicate if when using server-side processing the loading of data
		 * should be deferred until the second draw.
		 * Note that this parameter will be set by the initialisation routine. To
		 * set a default use {@link DataTable.defaults}.
		 *  @type boolean
		 *  @default false
		 */
		"bDeferLoading": false,
	
		/**
		 * Indicate if all required information has been read in
		 *  @type boolean
		 *  @default false
		 */
		"bInitialised": false,
	
		/**
		 * Information about open rows. Each object in the array has the parameters
		 * 'nTr' and 'nParent'
		 *  @type array
		 *  @default []
		 */
		"aoOpenRows": [],
	
		/**
		 * Dictate the positioning of DataTables' control elements - see
		 * {@link DataTable.model.oInit.sDom}.
		 * Note that this parameter will be set by the initialisation routine. To
		 * set a default use {@link DataTable.defaults}.
		 *  @type string
		 *  @default null
		 */
		"sDom": null,
	
		/**
		 * Search delay (in mS)
		 *  @type integer
		 *  @default null
		 */
		"searchDelay": null,
	
		/**
		 * Which type of pagination should be used.
		 * Note that this parameter will be set by the initialisation routine. To
		 * set a default use {@link DataTable.defaults}.
		 *  @type string
		 *  @default two_button
		 */
		"sPaginationType": "two_button",
	
		/**
		 * The state duration (for `stateSave`) in seconds.
		 * Note that this parameter will be set by the initialisation routine. To
		 * set a default use {@link DataTable.defaults}.
		 *  @type int
		 *  @default 0
		 */
		"iStateDuration": 0,
	
		/**
		 * Array of callback functions for state saving. Each array element is an
		 * object with the following parameters:
		 *   <ul>
		 *     <li>function:fn - function to call. Takes two parameters, oSettings
		 *       and the JSON string to save that has been thus far created. Returns
		 *       a JSON string to be inserted into a json object
		 *       (i.e. '"param": [ 0, 1, 2]')</li>
		 *     <li>string:sName - name of callback</li>
		 *   </ul>
		 *  @type array
		 *  @default []
		 */
		"aoStateSave": [],
	
		/**
		 * Array of callback functions for state loading. Each array element is an
		 * object with the following parameters:
		 *   <ul>
		 *     <li>function:fn - function to call. Takes two parameters, oSettings
		 *       and the object stored. May return false to cancel state loading</li>
		 *     <li>string:sName - name of callback</li>
		 *   </ul>
		 *  @type array
		 *  @default []
		 */
		"aoStateLoad": [],
	
		/**
		 * State that was saved. Useful for back reference
		 *  @type object
		 *  @default null
		 */
		"oSavedState": null,
	
		/**
		 * State that was loaded. Useful for back reference
		 *  @type object
		 *  @default null
		 */
		"oLoadedState": null,
	
		/**
		 * Source url for AJAX data for the table.
		 * Note that this parameter will be set by the initialisation routine. To
		 * set a default use {@link DataTable.defaults}.
		 *  @type string
		 *  @default null
		 */
		"sAjaxSource": null,
	
		/**
		 * Property from a given object from which to read the table data from. This
		 * can be an empty string (when not server-side processing), in which case
		 * it is  assumed an an array is given directly.
		 * Note that this parameter will be set by the initialisation routine. To
		 * set a default use {@link DataTable.defaults}.
		 *  @type string
		 */
		"sAjaxDataProp": null,
	
		/**
		 * The last jQuery XHR object that was used for server-side data gathering.
		 * This can be used for working with the XHR information in one of the
		 * callbacks
		 *  @type object
		 *  @default null
		 */
		"jqXHR": null,
	
		/**
		 * JSON returned from the server in the last Ajax request
		 *  @type object
		 *  @default undefined
		 */
		"json": undefined,
	
		/**
		 * Data submitted as part of the last Ajax request
		 *  @type object
		 *  @default undefined
		 */
		"oAjaxData": undefined,
	
		/**
		 * Function to get the server-side data.
		 * Note that this parameter will be set by the initialisation routine. To
		 * set a default use {@link DataTable.defaults}.
		 *  @type function
		 */
		"fnServerData": null,
	
		/**
		 * Functions which are called prior to sending an Ajax request so extra
		 * parameters can easily be sent to the server
		 *  @type array
		 *  @default []
		 */
		"aoServerParams": [],
	
		/**
		 * Send the XHR HTTP method - GET or POST (could be PUT or DELETE if
		 * required).
		 * Note that this parameter will be set by the initialisation routine. To
		 * set a default use {@link DataTable.defaults}.
		 *  @type string
		 */
		"sServerMethod": null,
	
		/**
		 * Format numbers for display.
		 * Note that this parameter will be set by the initialisation routine. To
		 * set a default use {@link DataTable.defaults}.
		 *  @type function
		 */
		"fnFormatNumber": null,
	
		/**
		 * List of options that can be used for the user selectable length menu.
		 * Note that this parameter will be set by the initialisation routine. To
		 * set a default use {@link DataTable.defaults}.
		 *  @type array
		 *  @default []
		 */
		"aLengthMenu": null,
	
		/**
		 * Counter for the draws that the table does. Also used as a tracker for
		 * server-side processing
		 *  @type int
		 *  @default 0
		 */
		"iDraw": 0,
	
		/**
		 * Indicate if a redraw is being done - useful for Ajax
		 *  @type boolean
		 *  @default false
		 */
		"bDrawing": false,
	
		/**
		 * Draw index (iDraw) of the last error when parsing the returned data
		 *  @type int
		 *  @default -1
		 */
		"iDrawError": -1,
	
		/**
		 * Paging display length
		 *  @type int
		 *  @default 10
		 */
		"_iDisplayLength": 10,
	
		/**
		 * Paging start point - aiDisplay index
		 *  @type int
		 *  @default 0
		 */
		"_iDisplayStart": 0,
	
		/**
		 * Server-side processing - number of records in the result set
		 * (i.e. before filtering), Use fnRecordsTotal rather than
		 * this property to get the value of the number of records, regardless of
		 * the server-side processing setting.
		 *  @type int
		 *  @default 0
		 *  @private
		 */
		"_iRecordsTotal": 0,
	
		/**
		 * Server-side processing - number of records in the current display set
		 * (i.e. after filtering). Use fnRecordsDisplay rather than
		 * this property to get the value of the number of records, regardless of
		 * the server-side processing setting.
		 *  @type boolean
		 *  @default 0
		 *  @private
		 */
		"_iRecordsDisplay": 0,
	
		/**
		 * The classes to use for the table
		 *  @type object
		 *  @default {}
		 */
		"oClasses": {},
	
		/**
		 * Flag attached to the settings object so you can check in the draw
		 * callback if filtering has been done in the draw. Deprecated in favour of
		 * events.
		 *  @type boolean
		 *  @default false
		 *  @deprecated
		 */
		"bFiltered": false,
	
		/**
		 * Flag attached to the settings object so you can check in the draw
		 * callback if sorting has been done in the draw. Deprecated in favour of
		 * events.
		 *  @type boolean
		 *  @default false
		 *  @deprecated
		 */
		"bSorted": false,
	
		/**
		 * Indicate that if multiple rows are in the header and there is more than
		 * one unique cell per column, if the top one (true) or bottom one (false)
		 * should be used for sorting / title by DataTables.
		 * Note that this parameter will be set by the initialisation routine. To
		 * set a default use {@link DataTable.defaults}.
		 *  @type boolean
		 */
		"bSortCellsTop": null,
	
		/**
		 * Initialisation object that is used for the table
		 *  @type object
		 *  @default null
		 */
		"oInit": null,
	
		/**
		 * Destroy callback functions - for plug-ins to attach themselves to the
		 * destroy so they can clean up markup and events.
		 *  @type array
		 *  @default []
		 */
		"aoDestroyCallback": [],
	
	
		/**
		 * Get the number of records in the current record set, before filtering
		 *  @type function
		 */
		"fnRecordsTotal": function ()
		{
			return _fnDataSource( this ) == 'ssp' ?
				this._iRecordsTotal * 1 :
				this.aiDisplayMaster.length;
		},
	
		/**
		 * Get the number of records in the current record set, after filtering
		 *  @type function
		 */
		"fnRecordsDisplay": function ()
		{
			return _fnDataSource( this ) == 'ssp' ?
				this._iRecordsDisplay * 1 :
				this.aiDisplay.length;
		},
	
		/**
		 * Get the display end point - aiDisplay index
		 *  @type function
		 */
		"fnDisplayEnd": function ()
		{
			var
				len      = this._iDisplayLength,
				start    = this._iDisplayStart,
				calc     = start + len,
				records  = this.aiDisplay.length,
				features = this.oFeatures,
				paginate = features.bPaginate;
	
			if ( features.bServerSide ) {
				return paginate === false || len === -1 ?
					start + records :
					Math.min( start+len, this._iRecordsDisplay );
			}
			else {
				return ! paginate || calc>records || len===-1 ?
					records :
					calc;
			}
		},
	
		/**
		 * The DataTables object for this table
		 *  @type object
		 *  @default null
		 */
		"oInstance": null,
	
		/**
		 * Unique identifier for each instance of the DataTables object. If there
		 * is an ID on the table node, then it takes that value, otherwise an
		 * incrementing internal counter is used.
		 *  @type string
		 *  @default null
		 */
		"sInstance": null,
	
		/**
		 * tabindex attribute value that is added to DataTables control elements, allowing
		 * keyboard navigation of the table and its controls.
		 */
		"iTabIndex": 0,
	
		/**
		 * DIV container for the footer scrolling table if scrolling
		 */
		"nScrollHead": null,
	
		/**
		 * DIV container for the footer scrolling table if scrolling
		 */
		"nScrollFoot": null,
	
		/**
		 * Last applied sort
		 *  @type array
		 *  @default []
		 */
		"aLastSort": [],
	
		/**
		 * Stored plug-in instances
		 *  @type object
		 *  @default {}
		 */
		"oPlugins": {},
	
		/**
		 * Function used to get a row's id from the row's data
		 *  @type function
		 *  @default null
		 */
		"rowIdFn": null,
	
		/**
		 * Data location where to store a row's id
		 *  @type string
		 *  @default null
		 */
		"rowId": null
	};
	
	/**
	 * Extension object for DataTables that is used to provide all extension
	 * options.
	 *
	 * Note that the `DataTable.ext` object is available through
	 * `jQuery.fn.dataTable.ext` where it may be accessed and manipulated. It is
	 * also aliased to `jQuery.fn.dataTableExt` for historic reasons.
	 *  @namespace
	 *  @extends DataTable.models.ext
	 */
	
	
	/**
	 * DataTables extensions
	 * 
	 * This namespace acts as a collection area for plug-ins that can be used to
	 * extend DataTables capabilities. Indeed many of the build in methods
	 * use this method to provide their own capabilities (sorting methods for
	 * example).
	 *
	 * Note that this namespace is aliased to `jQuery.fn.dataTableExt` for legacy
	 * reasons
	 *
	 *  @namespace
	 */
	DataTable.ext = _ext = {
		/**
		 * Buttons. For use with the Buttons extension for DataTables. This is
		 * defined here so other extensions can define buttons regardless of load
		 * order. It is _not_ used by DataTables core.
		 *
		 *  @type object
		 *  @default {}
		 */
		buttons: {},
	
	
		/**
		 * Element class names
		 *
		 *  @type object
		 *  @default {}
		 */
		classes: {},
	
	
		/**
		 * DataTables build type (expanded by the download builder)
		 *
		 *  @type string
		 */
		build:"bs5/dt-1.13.4/b-2.3.6/sb-1.4.2/sp-2.1.2",
	
	
		/**
		 * Error reporting.
		 * 
		 * How should DataTables report an error. Can take the value 'alert',
		 * 'throw', 'none' or a function.
		 *
		 *  @type string|function
		 *  @default alert
		 */
		errMode: "alert",
	
	
		/**
		 * Feature plug-ins.
		 * 
		 * This is an array of objects which describe the feature plug-ins that are
		 * available to DataTables. These feature plug-ins are then available for
		 * use through the `dom` initialisation option.
		 * 
		 * Each feature plug-in is described by an object which must have the
		 * following properties:
		 * 
		 * * `fnInit` - function that is used to initialise the plug-in,
		 * * `cFeature` - a character so the feature can be enabled by the `dom`
		 *   instillation option. This is case sensitive.
		 *
		 * The `fnInit` function has the following input parameters:
		 *
		 * 1. `{object}` DataTables settings object: see
		 *    {@link DataTable.models.oSettings}
		 *
		 * And the following return is expected:
		 * 
		 * * {node|null} The element which contains your feature. Note that the
		 *   return may also be void if your plug-in does not require to inject any
		 *   DOM elements into DataTables control (`dom`) - for example this might
		 *   be useful when developing a plug-in which allows table control via
		 *   keyboard entry
		 *
		 *  @type array
		 *
		 *  @example
		 *    $.fn.dataTable.ext.features.push( {
		 *      "fnInit": function( oSettings ) {
		 *        return new TableTools( { "oDTSettings": oSettings } );
		 *      },
		 *      "cFeature": "T"
		 *    } );
		 */
		feature: [],
	
	
		/**
		 * Row searching.
		 * 
		 * This method of searching is complimentary to the default type based
		 * searching, and a lot more comprehensive as it allows you complete control
		 * over the searching logic. Each element in this array is a function
		 * (parameters described below) that is called for every row in the table,
		 * and your logic decides if it should be included in the searching data set
		 * or not.
		 *
		 * Searching functions have the following input parameters:
		 *
		 * 1. `{object}` DataTables settings object: see
		 *    {@link DataTable.models.oSettings}
		 * 2. `{array|object}` Data for the row to be processed (same as the
		 *    original format that was passed in as the data source, or an array
		 *    from a DOM data source
		 * 3. `{int}` Row index ({@link DataTable.models.oSettings.aoData}), which
		 *    can be useful to retrieve the `TR` element if you need DOM interaction.
		 *
		 * And the following return is expected:
		 *
		 * * {boolean} Include the row in the searched result set (true) or not
		 *   (false)
		 *
		 * Note that as with the main search ability in DataTables, technically this
		 * is "filtering", since it is subtractive. However, for consistency in
		 * naming we call it searching here.
		 *
		 *  @type array
		 *  @default []
		 *
		 *  @example
		 *    // The following example shows custom search being applied to the
		 *    // fourth column (i.e. the data[3] index) based on two input values
		 *    // from the end-user, matching the data in a certain range.
		 *    $.fn.dataTable.ext.search.push(
		 *      function( settings, data, dataIndex ) {
		 *        var min = document.getElementById('min').value * 1;
		 *        var max = document.getElementById('max').value * 1;
		 *        var version = data[3] == "-" ? 0 : data[3]*1;
		 *
		 *        if ( min == "" && max == "" ) {
		 *          return true;
		 *        }
		 *        else if ( min == "" && version < max ) {
		 *          return true;
		 *        }
		 *        else if ( min < version && "" == max ) {
		 *          return true;
		 *        }
		 *        else if ( min < version && version < max ) {
		 *          return true;
		 *        }
		 *        return false;
		 *      }
		 *    );
		 */
		search: [],
	
	
		/**
		 * Selector extensions
		 *
		 * The `selector` option can be used to extend the options available for the
		 * selector modifier options (`selector-modifier` object data type) that
		 * each of the three built in selector types offer (row, column and cell +
		 * their plural counterparts). For example the Select extension uses this
		 * mechanism to provide an option to select only rows, columns and cells
		 * that have been marked as selected by the end user (`{selected: true}`),
		 * which can be used in conjunction with the existing built in selector
		 * options.
		 *
		 * Each property is an array to which functions can be pushed. The functions
		 * take three attributes:
		 *
		 * * Settings object for the host table
		 * * Options object (`selector-modifier` object type)
		 * * Array of selected item indexes
		 *
		 * The return is an array of the resulting item indexes after the custom
		 * selector has been applied.
		 *
		 *  @type object
		 */
		selector: {
			cell: [],
			column: [],
			row: []
		},
	
	
		/**
		 * Internal functions, exposed for used in plug-ins.
		 * 
		 * Please note that you should not need to use the internal methods for
		 * anything other than a plug-in (and even then, try to avoid if possible).
		 * The internal function may change between releases.
		 *
		 *  @type object
		 *  @default {}
		 */
		internal: {},
	
	
		/**
		 * Legacy configuration options. Enable and disable legacy options that
		 * are available in DataTables.
		 *
		 *  @type object
		 */
		legacy: {
			/**
			 * Enable / disable DataTables 1.9 compatible server-side processing
			 * requests
			 *
			 *  @type boolean
			 *  @default null
			 */
			ajax: null
		},
	
	
		/**
		 * Pagination plug-in methods.
		 * 
		 * Each entry in this object is a function and defines which buttons should
		 * be shown by the pagination rendering method that is used for the table:
		 * {@link DataTable.ext.renderer.pageButton}. The renderer addresses how the
		 * buttons are displayed in the document, while the functions here tell it
		 * what buttons to display. This is done by returning an array of button
		 * descriptions (what each button will do).
		 *
		 * Pagination types (the four built in options and any additional plug-in
		 * options defined here) can be used through the `paginationType`
		 * initialisation parameter.
		 *
		 * The functions defined take two parameters:
		 *
		 * 1. `{int} page` The current page index
		 * 2. `{int} pages` The number of pages in the table
		 *
		 * Each function is expected to return an array where each element of the
		 * array can be one of:
		 *
		 * * `first` - Jump to first page when activated
		 * * `last` - Jump to last page when activated
		 * * `previous` - Show previous page when activated
		 * * `next` - Show next page when activated
		 * * `{int}` - Show page of the index given
		 * * `{array}` - A nested array containing the above elements to add a
		 *   containing 'DIV' element (might be useful for styling).
		 *
		 * Note that DataTables v1.9- used this object slightly differently whereby
		 * an object with two functions would be defined for each plug-in. That
		 * ability is still supported by DataTables 1.10+ to provide backwards
		 * compatibility, but this option of use is now decremented and no longer
		 * documented in DataTables 1.10+.
		 *
		 *  @type object
		 *  @default {}
		 *
		 *  @example
		 *    // Show previous, next and current page buttons only
		 *    $.fn.dataTableExt.oPagination.current = function ( page, pages ) {
		 *      return [ 'previous', page, 'next' ];
		 *    };
		 */
		pager: {},
	
	
		renderer: {
			pageButton: {},
			header: {}
		},
	
	
		/**
		 * Ordering plug-ins - custom data source
		 * 
		 * The extension options for ordering of data available here is complimentary
		 * to the default type based ordering that DataTables typically uses. It
		 * allows much greater control over the the data that is being used to
		 * order a column, but is necessarily therefore more complex.
		 * 
		 * This type of ordering is useful if you want to do ordering based on data
		 * live from the DOM (for example the contents of an 'input' element) rather
		 * than just the static string that DataTables knows of.
		 * 
		 * The way these plug-ins work is that you create an array of the values you
		 * wish to be ordering for the column in question and then return that
		 * array. The data in the array much be in the index order of the rows in
		 * the table (not the currently ordering order!). Which order data gathering
		 * function is run here depends on the `dt-init columns.orderDataType`
		 * parameter that is used for the column (if any).
		 *
		 * The functions defined take two parameters:
		 *
		 * 1. `{object}` DataTables settings object: see
		 *    {@link DataTable.models.oSettings}
		 * 2. `{int}` Target column index
		 *
		 * Each function is expected to return an array:
		 *
		 * * `{array}` Data for the column to be ordering upon
		 *
		 *  @type array
		 *
		 *  @example
		 *    // Ordering using `input` node values
		 *    $.fn.dataTable.ext.order['dom-text'] = function  ( settings, col )
		 *    {
		 *      return this.api().column( col, {order:'index'} ).nodes().map( function ( td, i ) {
		 *        return $('input', td).val();
		 *      } );
		 *    }
		 */
		order: {},
	
	
		/**
		 * Type based plug-ins.
		 *
		 * Each column in DataTables has a type assigned to it, either by automatic
		 * detection or by direct assignment using the `type` option for the column.
		 * The type of a column will effect how it is ordering and search (plug-ins
		 * can also make use of the column type if required).
		 *
		 * @namespace
		 */
		type: {
			/**
			 * Type detection functions.
			 *
			 * The functions defined in this object are used to automatically detect
			 * a column's type, making initialisation of DataTables super easy, even
			 * when complex data is in the table.
			 *
			 * The functions defined take two parameters:
			 *
		     *  1. `{*}` Data from the column cell to be analysed
		     *  2. `{settings}` DataTables settings object. This can be used to
		     *     perform context specific type detection - for example detection
		     *     based on language settings such as using a comma for a decimal
		     *     place. Generally speaking the options from the settings will not
		     *     be required
			 *
			 * Each function is expected to return:
			 *
			 * * `{string|null}` Data type detected, or null if unknown (and thus
			 *   pass it on to the other type detection functions.
			 *
			 *  @type array
			 *
			 *  @example
			 *    // Currency type detection plug-in:
			 *    $.fn.dataTable.ext.type.detect.push(
			 *      function ( data, settings ) {
			 *        // Check the numeric part
			 *        if ( ! data.substring(1).match(/[0-9]/) ) {
			 *          return null;
			 *        }
			 *
			 *        // Check prefixed by currency
			 *        if ( data.charAt(0) == '$' || data.charAt(0) == '&pound;' ) {
			 *          return 'currency';
			 *        }
			 *        return null;
			 *      }
			 *    );
			 */
			detect: [],
	
	
			/**
			 * Type based search formatting.
			 *
			 * The type based searching functions can be used to pre-format the
			 * data to be search on. For example, it can be used to strip HTML
			 * tags or to de-format telephone numbers for numeric only searching.
			 *
			 * Note that is a search is not defined for a column of a given type,
			 * no search formatting will be performed.
			 * 
			 * Pre-processing of searching data plug-ins - When you assign the sType
			 * for a column (or have it automatically detected for you by DataTables
			 * or a type detection plug-in), you will typically be using this for
			 * custom sorting, but it can also be used to provide custom searching
			 * by allowing you to pre-processing the data and returning the data in
			 * the format that should be searched upon. This is done by adding
			 * functions this object with a parameter name which matches the sType
			 * for that target column. This is the corollary of <i>afnSortData</i>
			 * for searching data.
			 *
			 * The functions defined take a single parameter:
			 *
		     *  1. `{*}` Data from the column cell to be prepared for searching
			 *
			 * Each function is expected to return:
			 *
			 * * `{string|null}` Formatted string that will be used for the searching.
			 *
			 *  @type object
			 *  @default {}
			 *
			 *  @example
			 *    $.fn.dataTable.ext.type.search['title-numeric'] = function ( d ) {
			 *      return d.replace(/\n/g," ").replace( /<.*?>/g, "" );
			 *    }
			 */
			search: {},
	
	
			/**
			 * Type based ordering.
			 *
			 * The column type tells DataTables what ordering to apply to the table
			 * when a column is sorted upon. The order for each type that is defined,
			 * is defined by the functions available in this object.
			 *
			 * Each ordering option can be described by three properties added to
			 * this object:
			 *
			 * * `{type}-pre` - Pre-formatting function
			 * * `{type}-asc` - Ascending order function
			 * * `{type}-desc` - Descending order function
			 *
			 * All three can be used together, only `{type}-pre` or only
			 * `{type}-asc` and `{type}-desc` together. It is generally recommended
			 * that only `{type}-pre` is used, as this provides the optimal
			 * implementation in terms of speed, although the others are provided
			 * for compatibility with existing Javascript sort functions.
			 *
			 * `{type}-pre`: Functions defined take a single parameter:
			 *
		     *  1. `{*}` Data from the column cell to be prepared for ordering
			 *
			 * And return:
			 *
			 * * `{*}` Data to be sorted upon
			 *
			 * `{type}-asc` and `{type}-desc`: Functions are typical Javascript sort
			 * functions, taking two parameters:
			 *
		     *  1. `{*}` Data to compare to the second parameter
		     *  2. `{*}` Data to compare to the first parameter
			 *
			 * And returning:
			 *
			 * * `{*}` Ordering match: <0 if first parameter should be sorted lower
			 *   than the second parameter, ===0 if the two parameters are equal and
			 *   >0 if the first parameter should be sorted height than the second
			 *   parameter.
			 * 
			 *  @type object
			 *  @default {}
			 *
			 *  @example
			 *    // Numeric ordering of formatted numbers with a pre-formatter
			 *    $.extend( $.fn.dataTable.ext.type.order, {
			 *      "string-pre": function(x) {
			 *        a = (a === "-" || a === "") ? 0 : a.replace( /[^\d\-\.]/g, "" );
			 *        return parseFloat( a );
			 *      }
			 *    } );
			 *
			 *  @example
			 *    // Case-sensitive string ordering, with no pre-formatting method
			 *    $.extend( $.fn.dataTable.ext.order, {
			 *      "string-case-asc": function(x,y) {
			 *        return ((x < y) ? -1 : ((x > y) ? 1 : 0));
			 *      },
			 *      "string-case-desc": function(x,y) {
			 *        return ((x < y) ? 1 : ((x > y) ? -1 : 0));
			 *      }
			 *    } );
			 */
			order: {}
		},
	
		/**
		 * Unique DataTables instance counter
		 *
		 * @type int
		 * @private
		 */
		_unique: 0,
	
	
		//
		// Depreciated
		// The following properties are retained for backwards compatibility only.
		// The should not be used in new projects and will be removed in a future
		// version
		//
	
		/**
		 * Version check function.
		 *  @type function
		 *  @depreciated Since 1.10
		 */
		fnVersionCheck: DataTable.fnVersionCheck,
	
	
		/**
		 * Index for what 'this' index API functions should use
		 *  @type int
		 *  @deprecated Since v1.10
		 */
		iApiIndex: 0,
	
	
		/**
		 * jQuery UI class container
		 *  @type object
		 *  @deprecated Since v1.10
		 */
		oJUIClasses: {},
	
	
		/**
		 * Software version
		 *  @type string
		 *  @deprecated Since v1.10
		 */
		sVersion: DataTable.version
	};
	
	
	//
	// Backwards compatibility. Alias to pre 1.10 Hungarian notation counter parts
	//
	$.extend( _ext, {
		afnFiltering: _ext.search,
		aTypes:       _ext.type.detect,
		ofnSearch:    _ext.type.search,
		oSort:        _ext.type.order,
		afnSortData:  _ext.order,
		aoFeatures:   _ext.feature,
		oApi:         _ext.internal,
		oStdClasses:  _ext.classes,
		oPagination:  _ext.pager
	} );
	
	
	$.extend( DataTable.ext.classes, {
		"sTable": "dataTable",
		"sNoFooter": "no-footer",
	
		/* Paging buttons */
		"sPageButton": "paginate_button",
		"sPageButtonActive": "current",
		"sPageButtonDisabled": "disabled",
	
		/* Striping classes */
		"sStripeOdd": "odd",
		"sStripeEven": "even",
	
		/* Empty row */
		"sRowEmpty": "dataTables_empty",
	
		/* Features */
		"sWrapper": "dataTables_wrapper",
		"sFilter": "dataTables_filter",
		"sInfo": "dataTables_info",
		"sPaging": "dataTables_paginate paging_", /* Note that the type is postfixed */
		"sLength": "dataTables_length",
		"sProcessing": "dataTables_processing",
	
		/* Sorting */
		"sSortAsc": "sorting_asc",
		"sSortDesc": "sorting_desc",
		"sSortable": "sorting", /* Sortable in both directions */
		"sSortableAsc": "sorting_desc_disabled",
		"sSortableDesc": "sorting_asc_disabled",
		"sSortableNone": "sorting_disabled",
		"sSortColumn": "sorting_", /* Note that an int is postfixed for the sorting order */
	
		/* Filtering */
		"sFilterInput": "",
	
		/* Page length */
		"sLengthSelect": "",
	
		/* Scrolling */
		"sScrollWrapper": "dataTables_scroll",
		"sScrollHead": "dataTables_scrollHead",
		"sScrollHeadInner": "dataTables_scrollHeadInner",
		"sScrollBody": "dataTables_scrollBody",
		"sScrollFoot": "dataTables_scrollFoot",
		"sScrollFootInner": "dataTables_scrollFootInner",
	
		/* Misc */
		"sHeaderTH": "",
		"sFooterTH": "",
	
		// Deprecated
		"sSortJUIAsc": "",
		"sSortJUIDesc": "",
		"sSortJUI": "",
		"sSortJUIAscAllowed": "",
		"sSortJUIDescAllowed": "",
		"sSortJUIWrapper": "",
		"sSortIcon": "",
		"sJUIHeader": "",
		"sJUIFooter": ""
	} );
	
	
	var extPagination = DataTable.ext.pager;
	
	function _numbers ( page, pages ) {
		var
			numbers = [],
			buttons = extPagination.numbers_length,
			half = Math.floor( buttons / 2 ),
			i = 1;
	
		if ( pages <= buttons ) {
			numbers = _range( 0, pages );
		}
		else if ( page <= half ) {
			numbers = _range( 0, buttons-2 );
			numbers.push( 'ellipsis' );
			numbers.push( pages-1 );
		}
		else if ( page >= pages - 1 - half ) {
			numbers = _range( pages-(buttons-2), pages );
			numbers.splice( 0, 0, 'ellipsis' ); // no unshift in ie6
			numbers.splice( 0, 0, 0 );
		}
		else {
			numbers = _range( page-half+2, page+half-1 );
			numbers.push( 'ellipsis' );
			numbers.push( pages-1 );
			numbers.splice( 0, 0, 'ellipsis' );
			numbers.splice( 0, 0, 0 );
		}
	
		numbers.DT_el = 'span';
		return numbers;
	}
	
	
	$.extend( extPagination, {
		simple: function ( page, pages ) {
			return [ 'previous', 'next' ];
		},
	
		full: function ( page, pages ) {
			return [  'first', 'previous', 'next', 'last' ];
		},
	
		numbers: function ( page, pages ) {
			return [ _numbers(page, pages) ];
		},
	
		simple_numbers: function ( page, pages ) {
			return [ 'previous', _numbers(page, pages), 'next' ];
		},
	
		full_numbers: function ( page, pages ) {
			return [ 'first', 'previous', _numbers(page, pages), 'next', 'last' ];
		},
		
		first_last_numbers: function (page, pages) {
	 		return ['first', _numbers(page, pages), 'last'];
	 	},
	
		// For testing and plug-ins to use
		_numbers: _numbers,
	
		// Number of number buttons (including ellipsis) to show. _Must be odd!_
		numbers_length: 7
	} );
	
	
	$.extend( true, DataTable.ext.renderer, {
		pageButton: {
			_: function ( settings, host, idx, buttons, page, pages ) {
				var classes = settings.oClasses;
				var lang = settings.oLanguage.oPaginate;
				var aria = settings.oLanguage.oAria.paginate || {};
				var btnDisplay, btnClass;
	
				var attach = function( container, buttons ) {
					var i, ien, node, button, tabIndex;
					var disabledClass = classes.sPageButtonDisabled;
					var clickHandler = function ( e ) {
						_fnPageChange( settings, e.data.action, true );
					};
	
					for ( i=0, ien=buttons.length ; i<ien ; i++ ) {
						button = buttons[i];
	
						if ( Array.isArray( button ) ) {
							var inner = $( '<'+(button.DT_el || 'div')+'/>' )
								.appendTo( container );
							attach( inner, button );
						}
						else {
							btnDisplay = null;
							btnClass = button;
							tabIndex = settings.iTabIndex;
	
							switch ( button ) {
								case 'ellipsis':
									container.append('<span class="ellipsis">&#x2026;</span>');
									break;
	
								case 'first':
									btnDisplay = lang.sFirst;
	
									if ( page === 0 ) {
										tabIndex = -1;
										btnClass += ' ' + disabledClass;
									}
									break;
	
								case 'previous':
									btnDisplay = lang.sPrevious;
	
									if ( page === 0 ) {
										tabIndex = -1;
										btnClass += ' ' + disabledClass;
									}
									break;
	
								case 'next':
									btnDisplay = lang.sNext;
	
									if ( pages === 0 || page === pages-1 ) {
										tabIndex = -1;
										btnClass += ' ' + disabledClass;
									}
									break;
	
								case 'last':
									btnDisplay = lang.sLast;
	
									if ( pages === 0 || page === pages-1 ) {
										tabIndex = -1;
										btnClass += ' ' + disabledClass;
									}
									break;
	
								default:
									btnDisplay = settings.fnFormatNumber( button + 1 );
									btnClass = page === button ?
										classes.sPageButtonActive : '';
									break;
							}
	
							if ( btnDisplay !== null ) {
								var tag = settings.oInit.pagingTag || 'a';
								var disabled = btnClass.indexOf(disabledClass) !== -1;
			
	
								node = $('<'+tag+'>', {
										'class': classes.sPageButton+' '+btnClass,
										'aria-controls': settings.sTableId,
										'aria-disabled': disabled ? 'true' : null,
										'aria-label': aria[ button ],
										'aria-role': 'link',
										'aria-current': btnClass === classes.sPageButtonActive ? 'page' : null,
										'data-dt-idx': button,
										'tabindex': tabIndex,
										'id': idx === 0 && typeof button === 'string' ?
											settings.sTableId +'_'+ button :
											null
									} )
									.html( btnDisplay )
									.appendTo( container );
	
								_fnBindAction(
									node, {action: button}, clickHandler
								);
							}
						}
					}
				};
	
				// IE9 throws an 'unknown error' if document.activeElement is used
				// inside an iframe or frame. Try / catch the error. Not good for
				// accessibility, but neither are frames.
				var activeEl;
	
				try {
					// Because this approach is destroying and recreating the paging
					// elements, focus is lost on the select button which is bad for
					// accessibility. So we want to restore focus once the draw has
					// completed
					activeEl = $(host).find(document.activeElement).data('dt-idx');
				}
				catch (e) {}
	
				attach( $(host).empty(), buttons );
	
				if ( activeEl !== undefined ) {
					$(host).find( '[data-dt-idx='+activeEl+']' ).trigger('focus');
				}
			}
		}
	} );
	
	
	
	// Built in type detection. See model.ext.aTypes for information about
	// what is required from this methods.
	$.extend( DataTable.ext.type.detect, [
		// Plain numbers - first since V8 detects some plain numbers as dates
		// e.g. Date.parse('55') (but not all, e.g. Date.parse('22')...).
		function ( d, settings )
		{
			var decimal = settings.oLanguage.sDecimal;
			return _isNumber( d, decimal ) ? 'num'+decimal : null;
		},
	
		// Dates (only those recognised by the browser's Date.parse)
		function ( d, settings )
		{
			// V8 tries _very_ hard to make a string passed into `Date.parse()`
			// valid, so we need to use a regex to restrict date formats. Use a
			// plug-in for anything other than ISO8601 style strings
			if ( d && !(d instanceof Date) && ! _re_date.test(d) ) {
				return null;
			}
			var parsed = Date.parse(d);
			return (parsed !== null && !isNaN(parsed)) || _empty(d) ? 'date' : null;
		},
	
		// Formatted numbers
		function ( d, settings )
		{
			var decimal = settings.oLanguage.sDecimal;
			return _isNumber( d, decimal, true ) ? 'num-fmt'+decimal : null;
		},
	
		// HTML numeric
		function ( d, settings )
		{
			var decimal = settings.oLanguage.sDecimal;
			return _htmlNumeric( d, decimal ) ? 'html-num'+decimal : null;
		},
	
		// HTML numeric, formatted
		function ( d, settings )
		{
			var decimal = settings.oLanguage.sDecimal;
			return _htmlNumeric( d, decimal, true ) ? 'html-num-fmt'+decimal : null;
		},
	
		// HTML (this is strict checking - there must be html)
		function ( d, settings )
		{
			return _empty( d ) || (typeof d === 'string' && d.indexOf('<') !== -1) ?
				'html' : null;
		}
	] );
	
	
	
	// Filter formatting functions. See model.ext.ofnSearch for information about
	// what is required from these methods.
	// 
	// Note that additional search methods are added for the html numbers and
	// html formatted numbers by `_addNumericSort()` when we know what the decimal
	// place is
	
	
	$.extend( DataTable.ext.type.search, {
		html: function ( data ) {
			return _empty(data) ?
				data :
				typeof data === 'string' ?
					data
						.replace( _re_new_lines, " " )
						.replace( _re_html, "" ) :
					'';
		},
	
		string: function ( data ) {
			return _empty(data) ?
				data :
				typeof data === 'string' ?
					data.replace( _re_new_lines, " " ) :
					data;
		}
	} );
	
	
	
	var __numericReplace = function ( d, decimalPlace, re1, re2 ) {
		if ( d !== 0 && (!d || d === '-') ) {
			return -Infinity;
		}
		
		let type = typeof d;
	
		if (type === 'number' || type === 'bigint') {
			return d;
		}
	
		// If a decimal place other than `.` is used, it needs to be given to the
		// function so we can detect it and replace with a `.` which is the only
		// decimal place Javascript recognises - it is not locale aware.
		if ( decimalPlace ) {
			d = _numToDecimal( d, decimalPlace );
		}
	
		if ( d.replace ) {
			if ( re1 ) {
				d = d.replace( re1, '' );
			}
	
			if ( re2 ) {
				d = d.replace( re2, '' );
			}
		}
	
		return d * 1;
	};
	
	
	// Add the numeric 'deformatting' functions for sorting and search. This is done
	// in a function to provide an easy ability for the language options to add
	// additional methods if a non-period decimal place is used.
	function _addNumericSort ( decimalPlace ) {
		$.each(
			{
				// Plain numbers
				"num": function ( d ) {
					return __numericReplace( d, decimalPlace );
				},
	
				// Formatted numbers
				"num-fmt": function ( d ) {
					return __numericReplace( d, decimalPlace, _re_formatted_numeric );
				},
	
				// HTML numeric
				"html-num": function ( d ) {
					return __numericReplace( d, decimalPlace, _re_html );
				},
	
				// HTML numeric, formatted
				"html-num-fmt": function ( d ) {
					return __numericReplace( d, decimalPlace, _re_html, _re_formatted_numeric );
				}
			},
			function ( key, fn ) {
				// Add the ordering method
				_ext.type.order[ key+decimalPlace+'-pre' ] = fn;
	
				// For HTML types add a search formatter that will strip the HTML
				if ( key.match(/^html\-/) ) {
					_ext.type.search[ key+decimalPlace ] = _ext.type.search.html;
				}
			}
		);
	}
	
	
	// Default sort methods
	$.extend( _ext.type.order, {
		// Dates
		"date-pre": function ( d ) {
			var ts = Date.parse( d );
			return isNaN(ts) ? -Infinity : ts;
		},
	
		// html
		"html-pre": function ( a ) {
			return _empty(a) ?
				'' :
				a.replace ?
					a.replace( /<.*?>/g, "" ).toLowerCase() :
					a+'';
		},
	
		// string
		"string-pre": function ( a ) {
			// This is a little complex, but faster than always calling toString,
			// http://jsperf.com/tostring-v-check
			return _empty(a) ?
				'' :
				typeof a === 'string' ?
					a.toLowerCase() :
					! a.toString ?
						'' :
						a.toString();
		},
	
		// string-asc and -desc are retained only for compatibility with the old
		// sort methods
		"string-asc": function ( x, y ) {
			return ((x < y) ? -1 : ((x > y) ? 1 : 0));
		},
	
		"string-desc": function ( x, y ) {
			return ((x < y) ? 1 : ((x > y) ? -1 : 0));
		}
	} );
	
	
	// Numeric sorting types - order doesn't matter here
	_addNumericSort( '' );
	
	
	$.extend( true, DataTable.ext.renderer, {
		header: {
			_: function ( settings, cell, column, classes ) {
				// No additional mark-up required
				// Attach a sort listener to update on sort - note that using the
				// `DT` namespace will allow the event to be removed automatically
				// on destroy, while the `dt` namespaced event is the one we are
				// listening for
				$(settings.nTable).on( 'order.dt.DT', function ( e, ctx, sorting, columns ) {
					if ( settings !== ctx ) { // need to check this this is the host
						return;               // table, not a nested one
					}
	
					var colIdx = column.idx;
	
					cell
						.removeClass(
							classes.sSortAsc +' '+
							classes.sSortDesc
						)
						.addClass( columns[ colIdx ] == 'asc' ?
							classes.sSortAsc : columns[ colIdx ] == 'desc' ?
								classes.sSortDesc :
								column.sSortingClass
						);
				} );
			},
	
			jqueryui: function ( settings, cell, column, classes ) {
				$('<div/>')
					.addClass( classes.sSortJUIWrapper )
					.append( cell.contents() )
					.append( $('<span/>')
						.addClass( classes.sSortIcon+' '+column.sSortingClassJUI )
					)
					.appendTo( cell );
	
				// Attach a sort listener to update on sort
				$(settings.nTable).on( 'order.dt.DT', function ( e, ctx, sorting, columns ) {
					if ( settings !== ctx ) {
						return;
					}
	
					var colIdx = column.idx;
	
					cell
						.removeClass( classes.sSortAsc +" "+classes.sSortDesc )
						.addClass( columns[ colIdx ] == 'asc' ?
							classes.sSortAsc : columns[ colIdx ] == 'desc' ?
								classes.sSortDesc :
								column.sSortingClass
						);
	
					cell
						.find( 'span.'+classes.sSortIcon )
						.removeClass(
							classes.sSortJUIAsc +" "+
							classes.sSortJUIDesc +" "+
							classes.sSortJUI +" "+
							classes.sSortJUIAscAllowed +" "+
							classes.sSortJUIDescAllowed
						)
						.addClass( columns[ colIdx ] == 'asc' ?
							classes.sSortJUIAsc : columns[ colIdx ] == 'desc' ?
								classes.sSortJUIDesc :
								column.sSortingClassJUI
						);
				} );
			}
		}
	} );
	
	/*
	 * Public helper functions. These aren't used internally by DataTables, or
	 * called by any of the options passed into DataTables, but they can be used
	 * externally by developers working with DataTables. They are helper functions
	 * to make working with DataTables a little bit easier.
	 */
	
	var __htmlEscapeEntities = function ( d ) {
		if (Array.isArray(d)) {
			d = d.join(',');
		}
	
		return typeof d === 'string' ?
			d
				.replace(/&/g, '&amp;')
				.replace(/</g, '&lt;')
				.replace(/>/g, '&gt;')
				.replace(/"/g, '&quot;') :
			d;
	};
	
	// Common logic for moment, luxon or a date action
	function __mld( dt, momentFn, luxonFn, dateFn, arg1 ) {
		if (window.moment) {
			return dt[momentFn]( arg1 );
		}
		else if (window.luxon) {
			return dt[luxonFn]( arg1 );
		}
		
		return dateFn ? dt[dateFn]( arg1 ) : dt;
	}
	
	
	var __mlWarning = false;
	function __mldObj (d, format, locale) {
		var dt;
	
		if (window.moment) {
			dt = window.moment.utc( d, format, locale, true );
	
			if (! dt.isValid()) {
				return null;
			}
		}
		else if (window.luxon) {
			dt = format && typeof d === 'string'
				? window.luxon.DateTime.fromFormat( d, format )
				: window.luxon.DateTime.fromISO( d );
	
			if (! dt.isValid) {
				return null;
			}
	
			dt.setLocale(locale);
		}
		else if (! format) {
			// No format given, must be ISO
			dt = new Date(d);
		}
		else {
			if (! __mlWarning) {
				alert('DataTables warning: Formatted date without Moment.js or Luxon - https://datatables.net/tn/17');
			}
	
			__mlWarning = true;
		}
	
		return dt;
	}
	
	// Wrapper for date, datetime and time which all operate the same way with the exception of
	// the output string for auto locale support
	function __mlHelper (localeString) {
		return function ( from, to, locale, def ) {
			// Luxon and Moment support
			// Argument shifting
			if ( arguments.length === 0 ) {
				locale = 'en';
				to = null; // means toLocaleString
				from = null; // means iso8601
			}
			else if ( arguments.length === 1 ) {
				locale = 'en';
				to = from;
				from = null;
			}
			else if ( arguments.length === 2 ) {
				locale = to;
				to = from;
				from = null;
			}
	
			var typeName = 'datetime-' + to;
	
			// Add type detection and sorting specific to this date format - we need to be able to identify
			// date type columns as such, rather than as numbers in extensions. Hence the need for this.
			if (! DataTable.ext.type.order[typeName]) {
				// The renderer will give the value to type detect as the type!
				DataTable.ext.type.detect.unshift(function (d) {
					return d === typeName ? typeName : false;
				});
	
				// The renderer gives us Moment, Luxon or Date obects for the sorting, all of which have a
				// `valueOf` which gives milliseconds epoch
				DataTable.ext.type.order[typeName + '-asc'] = function (a, b) {
					var x = a.valueOf();
					var y = b.valueOf();
	
					return x === y
						? 0
						: x < y
							? -1
							: 1;
				}
	
				DataTable.ext.type.order[typeName + '-desc'] = function (a, b) {
					var x = a.valueOf();
					var y = b.valueOf();
	
					return x === y
						? 0
						: x > y
							? -1
							: 1;
				}
			}
		
			return function ( d, type ) {
				// Allow for a default value
				if (d === null || d === undefined) {
					if (def === '--now') {
						// We treat everything as UTC further down, so no changes are
						// made, as such need to get the local date / time as if it were
						// UTC
						var local = new Date();
						d = new Date( Date.UTC(
							local.getFullYear(), local.getMonth(), local.getDate(),
							local.getHours(), local.getMinutes(), local.getSeconds()
						) );
					}
					else {
						d = '';
					}
				}
	
				if (type === 'type') {
					// Typing uses the type name for fast matching
					return typeName;
				}
	
				if (d === '') {
					return type !== 'sort'
						? ''
						: __mldObj('0000-01-01 00:00:00', null, locale);
				}
	
				// Shortcut. If `from` and `to` are the same, we are using the renderer to
				// format for ordering, not display - its already in the display format.
				if ( to !== null && from === to && type !== 'sort' && type !== 'type' && ! (d instanceof Date) ) {
					return d;
				}
	
				var dt = __mldObj(d, from, locale);
	
				if (dt === null) {
					return d;
				}
	
				if (type === 'sort') {
					return dt;
				}
				
				var formatted = to === null
					? __mld(dt, 'toDate', 'toJSDate', '')[localeString]()
					: __mld(dt, 'format', 'toFormat', 'toISOString', to);
	
				// XSS protection
				return type === 'display' ?
					__htmlEscapeEntities( formatted ) :
					formatted;
			};
		}
	}
	
	// Based on locale, determine standard number formatting
	// Fallback for legacy browsers is US English
	var __thousands = ',';
	var __decimal = '.';
	
	if (Intl) {
		try {
			var num = new Intl.NumberFormat().formatToParts(100000.1);
		
			for (var i=0 ; i<num.length ; i++) {
				if (num[i].type === 'group') {
					__thousands = num[i].value;
				}
				else if (num[i].type === 'decimal') {
					__decimal = num[i].value;
				}
			}
		}
		catch (e) {
			// noop
		}
	}
	
	// Formatted date time detection - use by declaring the formats you are going to use
	DataTable.datetime = function ( format, locale ) {
		var typeName = 'datetime-detect-' + format;
	
		if (! locale) {
			locale = 'en';
		}
	
		if (! DataTable.ext.type.order[typeName]) {
			DataTable.ext.type.detect.unshift(function (d) {
				var dt = __mldObj(d, format, locale);
				return d === '' || dt ? typeName : false;
			});
	
			DataTable.ext.type.order[typeName + '-pre'] = function (d) {
				return __mldObj(d, format, locale) || 0;
			}
		}
	}
	
	/**
	 * Helpers for `columns.render`.
	 *
	 * The options defined here can be used with the `columns.render` initialisation
	 * option to provide a display renderer. The following functions are defined:
	 *
	 * * `number` - Will format numeric data (defined by `columns.data`) for
	 *   display, retaining the original unformatted data for sorting and filtering.
	 *   It takes 5 parameters:
	 *   * `string` - Thousands grouping separator
	 *   * `string` - Decimal point indicator
	 *   * `integer` - Number of decimal points to show
	 *   * `string` (optional) - Prefix.
	 *   * `string` (optional) - Postfix (/suffix).
	 * * `text` - Escape HTML to help prevent XSS attacks. It has no optional
	 *   parameters.
	 *
	 * @example
	 *   // Column definition using the number renderer
	 *   {
	 *     data: "salary",
	 *     render: $.fn.dataTable.render.number( '\'', '.', 0, '$' )
	 *   }
	 *
	 * @namespace
	 */
	DataTable.render = {
		date: __mlHelper('toLocaleDateString'),
		datetime: __mlHelper('toLocaleString'),
		time: __mlHelper('toLocaleTimeString'),
		number: function ( thousands, decimal, precision, prefix, postfix ) {
			// Auto locale detection
			if (thousands === null || thousands === undefined) {
				thousands = __thousands;
			}
	
			if (decimal === null || decimal === undefined) {
				decimal = __decimal;
			}
	
			return {
				display: function ( d ) {
					if ( typeof d !== 'number' && typeof d !== 'string' ) {
						return d;
					}
	
					if (d === '' || d === null) {
						return d;
					}
	
					var negative = d < 0 ? '-' : '';
					var flo = parseFloat( d );
	
					// If NaN then there isn't much formatting that we can do - just
					// return immediately, escaping any HTML (this was supposed to
					// be a number after all)
					if ( isNaN( flo ) ) {
						return __htmlEscapeEntities( d );
					}
	
					flo = flo.toFixed( precision );
					d = Math.abs( flo );
	
					var intPart = parseInt( d, 10 );
					var floatPart = precision ?
						decimal+(d - intPart).toFixed( precision ).substring( 2 ):
						'';
	
					// If zero, then can't have a negative prefix
					if (intPart === 0 && parseFloat(floatPart) === 0) {
						negative = '';
					}
	
					return negative + (prefix||'') +
						intPart.toString().replace(
							/\B(?=(\d{3})+(?!\d))/g, thousands
						) +
						floatPart +
						(postfix||'');
				}
			};
		},
	
		text: function () {
			return {
				display: __htmlEscapeEntities,
				filter: __htmlEscapeEntities
			};
		}
	};
	
	
	/*
	 * This is really a good bit rubbish this method of exposing the internal methods
	 * publicly... - To be fixed in 2.0 using methods on the prototype
	 */
	
	
	/**
	 * Create a wrapper function for exporting an internal functions to an external API.
	 *  @param {string} fn API function name
	 *  @returns {function} wrapped function
	 *  @memberof DataTable#internal
	 */
	function _fnExternApiFunc (fn)
	{
		return function() {
			var args = [_fnSettingsFromNode( this[DataTable.ext.iApiIndex] )].concat(
				Array.prototype.slice.call(arguments)
			);
			return DataTable.ext.internal[fn].apply( this, args );
		};
	}
	
	
	/**
	 * Reference to internal functions for use by plug-in developers. Note that
	 * these methods are references to internal functions and are considered to be
	 * private. If you use these methods, be aware that they are liable to change
	 * between versions.
	 *  @namespace
	 */
	$.extend( DataTable.ext.internal, {
		_fnExternApiFunc: _fnExternApiFunc,
		_fnBuildAjax: _fnBuildAjax,
		_fnAjaxUpdate: _fnAjaxUpdate,
		_fnAjaxParameters: _fnAjaxParameters,
		_fnAjaxUpdateDraw: _fnAjaxUpdateDraw,
		_fnAjaxDataSrc: _fnAjaxDataSrc,
		_fnAddColumn: _fnAddColumn,
		_fnColumnOptions: _fnColumnOptions,
		_fnAdjustColumnSizing: _fnAdjustColumnSizing,
		_fnVisibleToColumnIndex: _fnVisibleToColumnIndex,
		_fnColumnIndexToVisible: _fnColumnIndexToVisible,
		_fnVisbleColumns: _fnVisbleColumns,
		_fnGetColumns: _fnGetColumns,
		_fnColumnTypes: _fnColumnTypes,
		_fnApplyColumnDefs: _fnApplyColumnDefs,
		_fnHungarianMap: _fnHungarianMap,
		_fnCamelToHungarian: _fnCamelToHungarian,
		_fnLanguageCompat: _fnLanguageCompat,
		_fnBrowserDetect: _fnBrowserDetect,
		_fnAddData: _fnAddData,
		_fnAddTr: _fnAddTr,
		_fnNodeToDataIndex: _fnNodeToDataIndex,
		_fnNodeToColumnIndex: _fnNodeToColumnIndex,
		_fnGetCellData: _fnGetCellData,
		_fnSetCellData: _fnSetCellData,
		_fnSplitObjNotation: _fnSplitObjNotation,
		_fnGetObjectDataFn: _fnGetObjectDataFn,
		_fnSetObjectDataFn: _fnSetObjectDataFn,
		_fnGetDataMaster: _fnGetDataMaster,
		_fnClearTable: _fnClearTable,
		_fnDeleteIndex: _fnDeleteIndex,
		_fnInvalidate: _fnInvalidate,
		_fnGetRowElements: _fnGetRowElements,
		_fnCreateTr: _fnCreateTr,
		_fnBuildHead: _fnBuildHead,
		_fnDrawHead: _fnDrawHead,
		_fnDraw: _fnDraw,
		_fnReDraw: _fnReDraw,
		_fnAddOptionsHtml: _fnAddOptionsHtml,
		_fnDetectHeader: _fnDetectHeader,
		_fnGetUniqueThs: _fnGetUniqueThs,
		_fnFeatureHtmlFilter: _fnFeatureHtmlFilter,
		_fnFilterComplete: _fnFilterComplete,
		_fnFilterCustom: _fnFilterCustom,
		_fnFilterColumn: _fnFilterColumn,
		_fnFilter: _fnFilter,
		_fnFilterCreateSearch: _fnFilterCreateSearch,
		_fnEscapeRegex: _fnEscapeRegex,
		_fnFilterData: _fnFilterData,
		_fnFeatureHtmlInfo: _fnFeatureHtmlInfo,
		_fnUpdateInfo: _fnUpdateInfo,
		_fnInfoMacros: _fnInfoMacros,
		_fnInitialise: _fnInitialise,
		_fnInitComplete: _fnInitComplete,
		_fnLengthChange: _fnLengthChange,
		_fnFeatureHtmlLength: _fnFeatureHtmlLength,
		_fnFeatureHtmlPaginate: _fnFeatureHtmlPaginate,
		_fnPageChange: _fnPageChange,
		_fnFeatureHtmlProcessing: _fnFeatureHtmlProcessing,
		_fnProcessingDisplay: _fnProcessingDisplay,
		_fnFeatureHtmlTable: _fnFeatureHtmlTable,
		_fnScrollDraw: _fnScrollDraw,
		_fnApplyToChildren: _fnApplyToChildren,
		_fnCalculateColumnWidths: _fnCalculateColumnWidths,
		_fnThrottle: _fnThrottle,
		_fnConvertToWidth: _fnConvertToWidth,
		_fnGetWidestNode: _fnGetWidestNode,
		_fnGetMaxLenString: _fnGetMaxLenString,
		_fnStringToCss: _fnStringToCss,
		_fnSortFlatten: _fnSortFlatten,
		_fnSort: _fnSort,
		_fnSortAria: _fnSortAria,
		_fnSortListener: _fnSortListener,
		_fnSortAttachListener: _fnSortAttachListener,
		_fnSortingClasses: _fnSortingClasses,
		_fnSortData: _fnSortData,
		_fnSaveState: _fnSaveState,
		_fnLoadState: _fnLoadState,
		_fnImplementState: _fnImplementState,
		_fnSettingsFromNode: _fnSettingsFromNode,
		_fnLog: _fnLog,
		_fnMap: _fnMap,
		_fnBindAction: _fnBindAction,
		_fnCallbackReg: _fnCallbackReg,
		_fnCallbackFire: _fnCallbackFire,
		_fnLengthOverflow: _fnLengthOverflow,
		_fnRenderer: _fnRenderer,
		_fnDataSource: _fnDataSource,
		_fnRowAttributes: _fnRowAttributes,
		_fnExtend: _fnExtend,
		_fnCalculateEnd: function () {} // Used by a lot of plug-ins, but redundant
		                                // in 1.10, so this dead-end function is
		                                // added to prevent errors
	} );
	
	
	// jQuery access
	$.fn.dataTable = DataTable;
	
	// Provide access to the host jQuery object (circular reference)
	DataTable.$ = $;
	
	// Legacy aliases
	$.fn.dataTableSettings = DataTable.settings;
	$.fn.dataTableExt = DataTable.ext;
	
	// With a capital `D` we return a DataTables API instance rather than a
	// jQuery object
	$.fn.DataTable = function ( opts ) {
		return $(this).dataTable( opts ).api();
	};
	
	// All properties that are available to $.fn.dataTable should also be
	// available on $.fn.DataTable
	$.each( DataTable, function ( prop, val ) {
		$.fn.DataTable[ prop ] = val;
	} );

	return DataTable;
}));


/*! DataTables Bootstrap 5 integration
 * 2020 SpryMedia Ltd - datatables.net/license
 */

(function( factory ){
	if ( typeof define === 'function' && define.amd ) {
		// AMD
		define( ['jquery', 'datatables.net'], function ( $ ) {
			return factory( $, window, document );
		} );
	}
	else if ( typeof exports === 'object' ) {
		// CommonJS
		var jq = require('jquery');
		var cjsRequires = function (root, $) {
			if ( ! $.fn.dataTable ) {
				require('datatables.net')(root, $);
			}
		};

		if (typeof window !== 'undefined') {
			module.exports = function (root, $) {
				if ( ! root ) {
					// CommonJS environments without a window global must pass a
					// root. This will give an error otherwise
					root = window;
				}

				if ( ! $ ) {
					$ = jq( root );
				}

				cjsRequires( root, $ );
				return factory( $, root, root.document );
			};
		}
		else {
			cjsRequires( window, jq );
			module.exports = factory( jq, window, window.document );
		}
	}
	else {
		// Browser
		factory( jQuery, window, document );
	}
}(function( $, window, document, undefined ) {
'use strict';
var DataTable = $.fn.dataTable;



/**
 * DataTables integration for Bootstrap 5. This requires Bootstrap 5 and
 * DataTables 1.10 or newer.
 *
 * This file sets the defaults and adds options to DataTables to style its
 * controls using Bootstrap. See http://datatables.net/manual/styling/bootstrap
 * for further information.
 */

/* Set the defaults for DataTables initialisation */
$.extend( true, DataTable.defaults, {
	dom:
		"<'row'<'col-sm-12 col-md-6'l><'col-sm-12 col-md-6'f>>" +
		"<'row dt-row'<'col-sm-12'tr>>" +
		"<'row'<'col-sm-12 col-md-5'i><'col-sm-12 col-md-7'p>>",
	renderer: 'bootstrap'
} );


/* Default class modification */
$.extend( DataTable.ext.classes, {
	sWrapper:      "dataTables_wrapper dt-bootstrap5",
	sFilterInput:  "form-control form-control-sm",
	sLengthSelect: "form-select form-select-sm",
	sProcessing:   "dataTables_processing card",
	sPageButton:   "paginate_button page-item"
} );


/* Bootstrap paging button renderer */
DataTable.ext.renderer.pageButton.bootstrap = function ( settings, host, idx, buttons, page, pages ) {
	var api     = new DataTable.Api( settings );
	var classes = settings.oClasses;
	var lang    = settings.oLanguage.oPaginate;
	var aria = settings.oLanguage.oAria.paginate || {};
	var btnDisplay, btnClass;

	var attach = function( container, buttons ) {
		var i, ien, node, button;
		var clickHandler = function ( e ) {
			e.preventDefault();
			if ( !$(e.currentTarget).hasClass('disabled') && api.page() != e.data.action ) {
				api.page( e.data.action ).draw( 'page' );
			}
		};

		for ( i=0, ien=buttons.length ; i<ien ; i++ ) {
			button = buttons[i];

			if ( Array.isArray( button ) ) {
				attach( container, button );
			}
			else {
				btnDisplay = '';
				btnClass = '';

				switch ( button ) {
					case 'ellipsis':
						btnDisplay = '&#x2026;';
						btnClass = 'disabled';
						break;

					case 'first':
						btnDisplay = lang.sFirst;
						btnClass = button + (page > 0 ?
							'' : ' disabled');
						break;

					case 'previous':
						btnDisplay = lang.sPrevious;
						btnClass = button + (page > 0 ?
							'' : ' disabled');
						break;

					case 'next':
						btnDisplay = lang.sNext;
						btnClass = button + (page < pages-1 ?
							'' : ' disabled');
						break;

					case 'last':
						btnDisplay = lang.sLast;
						btnClass = button + (page < pages-1 ?
							'' : ' disabled');
						break;

					default:
						btnDisplay = button + 1;
						btnClass = page === button ?
							'active' : '';
						break;
				}

				if ( btnDisplay ) {
					var disabled = btnClass.indexOf('disabled') !== -1;

					node = $('<li>', {
							'class': classes.sPageButton+' '+btnClass,
							'id': idx === 0 && typeof button === 'string' ?
								settings.sTableId +'_'+ button :
								null
						} )
						.append( $('<a>', {
								'href': disabled ? null : '#',
								'aria-controls': settings.sTableId,
								'aria-disabled': disabled ? 'true' : null,
								'aria-label': aria[ button ],
								'aria-role': 'link',
								'aria-current': btnClass === 'active' ? 'page' : null,
								'data-dt-idx': button,
								'tabindex': settings.iTabIndex,
								'class': 'page-link'
							} )
							.html( btnDisplay )
						)
						.appendTo( container );

					settings.oApi._fnBindAction(
						node, {action: button}, clickHandler
					);
				}
			}
		}
	};

	var hostEl = $(host);
	// IE9 throws an 'unknown error' if document.activeElement is used
	// inside an iframe or frame. 
	var activeEl;

	try {
		// Because this approach is destroying and recreating the paging
		// elements, focus is lost on the select button which is bad for
		// accessibility. So we want to restore focus once the draw has
		// completed
		activeEl = hostEl.find(document.activeElement).data('dt-idx');
	}
	catch (e) {}

	var paginationEl = hostEl.children('ul.pagination');

	if (paginationEl.length) {
		paginationEl.empty();
	}
	else {
		paginationEl = hostEl.html('<ul/>').children('ul').addClass('pagination');
	}

	attach(
		paginationEl,
		buttons
	);

	if ( activeEl !== undefined ) {
		hostEl.find('[data-dt-idx='+activeEl+']').trigger('focus');
	}
};


return DataTable;
}));


/*! Buttons for DataTables 2.3.6
 * ©2016-2023 SpryMedia Ltd - datatables.net/license
 */

(function( factory ){
	if ( typeof define === 'function' && define.amd ) {
		// AMD
		define( ['jquery', 'datatables.net'], function ( $ ) {
			return factory( $, window, document );
		} );
	}
	else if ( typeof exports === 'object' ) {
		// CommonJS
		var jq = require('jquery');
		var cjsRequires = function (root, $) {
			if ( ! $.fn.dataTable ) {
				require('datatables.net')(root, $);
			}
		};

		if (typeof window !== 'undefined') {
			module.exports = function (root, $) {
				if ( ! root ) {
					// CommonJS environments without a window global must pass a
					// root. This will give an error otherwise
					root = window;
				}

				if ( ! $ ) {
					$ = jq( root );
				}

				cjsRequires( root, $ );
				return factory( $, root, root.document );
			};
		}
		else {
			cjsRequires( window, jq );
			module.exports = factory( jq, window, window.document );
		}
	}
	else {
		// Browser
		factory( jQuery, window, document );
	}
}(function( $, window, document, undefined ) {
'use strict';
var DataTable = $.fn.dataTable;



// Used for namespacing events added to the document by each instance, so they
// can be removed on destroy
var _instCounter = 0;

// Button namespacing counter for namespacing events on individual buttons
var _buttonCounter = 0;

var _dtButtons = DataTable.ext.buttons;

// Allow for jQuery slim
function _fadeIn(el, duration, fn) {
	if ($.fn.animate) {
		el
			.stop()
			.fadeIn( duration, fn );

	}
	else {
		el.css('display', 'block');

		if (fn) {
			fn.call(el);
		}
	}
}

function _fadeOut(el, duration, fn) {
	if ($.fn.animate) {
		el
			.stop()
			.fadeOut( duration, fn );
	}
	else {
		el.css('display', 'none');
		
		if (fn) {
			fn.call(el);
		}
	}
}

/**
 * [Buttons description]
 * @param {[type]}
 * @param {[type]}
 */
var Buttons = function( dt, config )
{
	// If not created with a `new` keyword then we return a wrapper function that
	// will take the settings object for a DT. This allows easy use of new instances
	// with the `layout` option - e.g. `topLeft: $.fn.dataTable.Buttons( ... )`.
	if ( !(this instanceof Buttons) ) {
		return function (settings) {
			return new Buttons( settings, dt ).container();
		};
	}

	// If there is no config set it to an empty object
	if ( typeof( config ) === 'undefined' ) {
		config = {};	
	}
	
	// Allow a boolean true for defaults
	if ( config === true ) {
		config = {};
	}

	// For easy configuration of buttons an array can be given
	if ( Array.isArray( config ) ) {
		config = { buttons: config };
	}

	this.c = $.extend( true, {}, Buttons.defaults, config );

	// Don't want a deep copy for the buttons
	if ( config.buttons ) {
		this.c.buttons = config.buttons;
	}

	this.s = {
		dt: new DataTable.Api( dt ),
		buttons: [],
		listenKeys: '',
		namespace: 'dtb'+(_instCounter++)
	};

	this.dom = {
		container: $('<'+this.c.dom.container.tag+'/>')
			.addClass( this.c.dom.container.className )
	};

	this._constructor();
};


$.extend( Buttons.prototype, {
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Public methods
	 */

	/**
	 * Get the action of a button
	 * @param  {int|string} Button index
	 * @return {function}
	 *//**
	 * Set the action of a button
	 * @param  {node} node Button element
	 * @param  {function} action Function to set
	 * @return {Buttons} Self for chaining
	 */
	action: function ( node, action )
	{
		var button = this._nodeToButton( node );

		if ( action === undefined ) {
			return button.conf.action;
		}

		button.conf.action = action;

		return this;
	},

	/**
	 * Add an active class to the button to make to look active or get current
	 * active state.
	 * @param  {node} node Button element
	 * @param  {boolean} [flag] Enable / disable flag
	 * @return {Buttons} Self for chaining or boolean for getter
	 */
	active: function ( node, flag ) {
		var button = this._nodeToButton( node );
		var klass = this.c.dom.button.active;
		var jqNode = $(button.node);

		if ( flag === undefined ) {
			return jqNode.hasClass( klass );
		}

		jqNode.toggleClass( klass, flag === undefined ? true : flag );

		return this;
	},

	/**
	 * Add a new button
	 * @param {object} config Button configuration object, base string name or function
	 * @param {int|string} [idx] Button index for where to insert the button
	 * @param {boolean} [draw=true] Trigger a draw. Set a false when adding
	 *   lots of buttons, until the last button.
	 * @return {Buttons} Self for chaining
	 */
	add: function ( config, idx, draw )
	{
		var buttons = this.s.buttons;

		if ( typeof idx === 'string' ) {
			var split = idx.split('-');
			var base = this.s;

			for ( var i=0, ien=split.length-1 ; i<ien ; i++ ) {
				base = base.buttons[ split[i]*1 ];
			}

			buttons = base.buttons;
			idx = split[ split.length-1 ]*1;
		}

		this._expandButton(
			buttons,
			config,
			config !== undefined ? config.split : undefined,
			(config === undefined || config.split === undefined || config.split.length === 0) && base !== undefined,
			false,
			idx
		);

		if (draw === undefined || draw === true) {
			this._draw();
		}
	
		return this;
	},

	/**
	 * Clear buttons from a collection and then insert new buttons
	 */
	collectionRebuild: function ( node, newButtons )
	{
		var button = this._nodeToButton( node );
		
		if(newButtons !== undefined) {
			var i;
			// Need to reverse the array
			for (i=button.buttons.length-1; i>=0; i--) {
				this.remove(button.buttons[i].node);
			}

			// If the collection has prefix and / or postfix buttons we need to add them in
			if (button.conf.prefixButtons) {
				newButtons.unshift.apply(newButtons, button.conf.prefixButtons);
			}

			if (button.conf.postfixButtons) {
				newButtons.push.apply(newButtons, button.conf.postfixButtons);
			}

			for (i=0; i<newButtons.length; i++) {
				var newBtn = newButtons[i];

				this._expandButton(
					button.buttons,
					newBtn,
					newBtn !== undefined && newBtn.config !== undefined && newBtn.config.split !== undefined,
					true,
					newBtn.parentConf !== undefined && newBtn.parentConf.split !== undefined,
					null,
					newBtn.parentConf
				);
			}
		}

		this._draw(button.collection, button.buttons);
	},

	/**
	 * Get the container node for the buttons
	 * @return {jQuery} Buttons node
	 */
	container: function ()
	{
		return this.dom.container;
	},

	/**
	 * Disable a button
	 * @param  {node} node Button node
	 * @return {Buttons} Self for chaining
	 */
	disable: function ( node ) {
		var button = this._nodeToButton( node );

		$(button.node)
			.addClass( this.c.dom.button.disabled )
			.prop('disabled', true);

		return this;
	},

	/**
	 * Destroy the instance, cleaning up event handlers and removing DOM
	 * elements
	 * @return {Buttons} Self for chaining
	 */
	destroy: function ()
	{
		// Key event listener
		$('body').off( 'keyup.'+this.s.namespace );

		// Individual button destroy (so they can remove their own events if
		// needed). Take a copy as the array is modified by `remove`
		var buttons = this.s.buttons.slice();
		var i, ien;
		
		for ( i=0, ien=buttons.length ; i<ien ; i++ ) {
			this.remove( buttons[i].node );
		}

		// Container
		this.dom.container.remove();

		// Remove from the settings object collection
		var buttonInsts = this.s.dt.settings()[0];

		for ( i=0, ien=buttonInsts.length ; i<ien ; i++ ) {
			if ( buttonInsts.inst === this ) {
				buttonInsts.splice( i, 1 );
				break;
			}
		}

		return this;
	},

	/**
	 * Enable / disable a button
	 * @param  {node} node Button node
	 * @param  {boolean} [flag=true] Enable / disable flag
	 * @return {Buttons} Self for chaining
	 */
	enable: function ( node, flag )
	{
		if ( flag === false ) {
			return this.disable( node );
		}

		var button = this._nodeToButton( node );
		$(button.node)
			.removeClass( this.c.dom.button.disabled )
			.prop('disabled', false);

		return this;
	},

	/**
	 * Get a button's index
	 * 
	 * This is internally recursive
	 * @param {element} node Button to get the index of
	 * @return {string} Button index
	 */
	index: function ( node, nested, buttons )
	{
		if ( ! nested ) {
			nested = '';
			buttons = this.s.buttons;
		}

		for ( var i=0, ien=buttons.length ; i<ien ; i++ ) {
			var inner = buttons[i].buttons;

			if (buttons[i].node === node) {
				return nested + i;
			}

			if ( inner && inner.length ) {
				var match = this.index(node, i + '-', inner);

				if (match !== null) {
					return match;
				}
			}
		}

		return null;
	},


	/**
	 * Get the instance name for the button set selector
	 * @return {string} Instance name
	 */
	name: function ()
	{
		return this.c.name;
	},

	/**
	 * Get a button's node of the buttons container if no button is given
	 * @param  {node} [node] Button node
	 * @return {jQuery} Button element, or container
	 */
	node: function ( node )
	{
		if ( ! node ) {
			return this.dom.container;
		}

		var button = this._nodeToButton( node );
		return $(button.node);
	},

	/**
	 * Set / get a processing class on the selected button
	 * @param {element} node Triggering button node
	 * @param  {boolean} flag true to add, false to remove, undefined to get
	 * @return {boolean|Buttons} Getter value or this if a setter.
	 */
	processing: function ( node, flag )
	{
		var dt = this.s.dt;
		var button = this._nodeToButton( node );

		if ( flag === undefined ) {
			return $(button.node).hasClass( 'processing' );
		}

		$(button.node).toggleClass( 'processing', flag );

		$(dt.table().node()).triggerHandler( 'buttons-processing.dt', [
			flag, dt.button( node ), dt, $(node), button.conf
		] );

		return this;
	},

	/**
	 * Remove a button.
	 * @param  {node} node Button node
	 * @return {Buttons} Self for chaining
	 */
	remove: function ( node )
	{
		var button = this._nodeToButton( node );
		var host = this._nodeToHost( node );
		var dt = this.s.dt;

		// Remove any child buttons first
		if ( button.buttons.length ) {
			for ( var i=button.buttons.length-1 ; i>=0 ; i-- ) {
				this.remove( button.buttons[i].node );
			}
		}

		button.conf.destroying = true;

		// Allow the button to remove event handlers, etc
		if ( button.conf.destroy ) {
			button.conf.destroy.call( dt.button(node), dt, $(node), button.conf );
		}

		this._removeKey( button.conf );

		$(button.node).remove();

		var idx = $.inArray( button, host );
		host.splice( idx, 1 );

		return this;
	},

	/**
	 * Get the text for a button
	 * @param  {int|string} node Button index
	 * @return {string} Button text
	 *//**
	 * Set the text for a button
	 * @param  {int|string|function} node Button index
	 * @param  {string} label Text
	 * @return {Buttons} Self for chaining
	 */
	text: function ( node, label )
	{
		var button = this._nodeToButton( node );
		var buttonLiner = this.c.dom.collection.buttonLiner;
		var linerTag = button.inCollection && buttonLiner && buttonLiner.tag ?
			buttonLiner.tag :
			this.c.dom.buttonLiner.tag;
		var dt = this.s.dt;
		var jqNode = $(button.node);
		var text = function ( opt ) {
			return typeof opt === 'function' ?
				opt( dt, jqNode, button.conf ) :
				opt;
		};

		if ( label === undefined ) {
			return text( button.conf.text );
		}

		button.conf.text = label;

		if ( linerTag ) {
			jqNode
				.children( linerTag )
				.eq(0)
				.filter(':not(.dt-down-arrow)')
				.html( text(label) );
		}
		else {
			jqNode.html( text(label) );
		}

		return this;
	},


	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Constructor
	 */

	/**
	 * Buttons constructor
	 * @private
	 */
	_constructor: function ()
	{
		var that = this;
		var dt = this.s.dt;
		var dtSettings = dt.settings()[0];
		var buttons =  this.c.buttons;

		if ( ! dtSettings._buttons ) {
			dtSettings._buttons = [];
		}

		dtSettings._buttons.push( {
			inst: this,
			name: this.c.name
		} );

		for ( var i=0, ien=buttons.length ; i<ien ; i++ ) {
			this.add( buttons[i] );
		}

		dt.on( 'destroy', function ( e, settings ) {
			if ( settings === dtSettings ) {
				that.destroy();
			}
		} );

		// Global key event binding to listen for button keys
		$('body').on( 'keyup.'+this.s.namespace, function ( e ) {
			if ( ! document.activeElement || document.activeElement === document.body ) {
				// SUse a string of characters for fast lookup of if we need to
				// handle this
				var character = String.fromCharCode(e.keyCode).toLowerCase();

				if ( that.s.listenKeys.toLowerCase().indexOf( character ) !== -1 ) {
					that._keypress( character, e );
				}
			}
		} );
	},


	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Private methods
	 */

	/**
	 * Add a new button to the key press listener
	 * @param {object} conf Resolved button configuration object
	 * @private
	 */
	_addKey: function ( conf )
	{
		if ( conf.key ) {
			this.s.listenKeys += $.isPlainObject( conf.key ) ?
				conf.key.key :
				conf.key;
		}
	},

	/**
	 * Insert the buttons into the container. Call without parameters!
	 * @param  {node} [container] Recursive only - Insert point
	 * @param  {array} [buttons] Recursive only - Buttons array
	 * @private
	 */
	_draw: function ( container, buttons )
	{
		if ( ! container ) {
			container = this.dom.container;
			buttons = this.s.buttons;
		}

		container.children().detach();

		for ( var i=0, ien=buttons.length ; i<ien ; i++ ) {
			container.append( buttons[i].inserter );
			container.append( ' ' );

			if ( buttons[i].buttons && buttons[i].buttons.length ) {
				this._draw( buttons[i].collection, buttons[i].buttons );
			}
		}
	},

	/**
	 * Create buttons from an array of buttons
	 * @param  {array} attachTo Buttons array to attach to
	 * @param  {object} button Button definition
	 * @param  {boolean} inCollection true if the button is in a collection
	 * @private
	 */
	_expandButton: function ( attachTo, button, split, inCollection, inSplit, attachPoint, parentConf )
	{
		var dt = this.s.dt;
		var buttonCounter = 0;
		var isSplit = false;
		var buttons = ! Array.isArray( button ) ?
			[ button ] :
			button;
		
		if(button === undefined ) {
			buttons = !Array.isArray(split) ?
				[ split ] :
				split;
		}

		if (button !== undefined && button.split !== undefined) {
			isSplit = true;
		}
			
		for ( var i=0, ien=buttons.length ; i<ien ; i++ ) {
			var conf = this._resolveExtends( buttons[i] );

			if ( ! conf ) {
				continue;
			}

			if( conf.config !== undefined && conf.config.split) {
				isSplit = true;
			}
			else {
				isSplit = false;
			}
			
			// If the configuration is an array, then expand the buttons at this
			// point
			if ( Array.isArray( conf ) ) {
				this._expandButton( attachTo, conf, built !== undefined && built.conf !== undefined ? built.conf.split : undefined, inCollection, parentConf !== undefined && parentConf.split !== undefined, attachPoint, parentConf );
				continue;
			}

			var built = this._buildButton( conf, inCollection, conf.split !== undefined || (conf.config !== undefined && conf.config.split !== undefined), inSplit );
			if ( ! built ) {
				continue;
			}

			if ( attachPoint !== undefined && attachPoint !== null ) {
				attachTo.splice( attachPoint, 0, built );
				attachPoint++;
			}
			else {
				attachTo.push( built );
			}

			
			if ( built.conf.buttons || built.conf.split ) {
				built.collection = $('<'+(isSplit ? this.c.dom.splitCollection.tag : this.c.dom.collection.tag)+'/>');

				built.conf._collection = built.collection;

				if(built.conf.split) {
					for(var j = 0; j < built.conf.split.length; j++) {
						if(typeof built.conf.split[j] === "object") {
							built.conf.split[j].parent = parentConf;
							if(built.conf.split[j].collectionLayout === undefined) {
								built.conf.split[j].collectionLayout = built.conf.collectionLayout;
							}
							if(built.conf.split[j].dropup === undefined) {
								built.conf.split[j].dropup = built.conf.dropup;
							}
							if(built.conf.split[j].fade === undefined) {
								built.conf.split[j].fade = built.conf.fade;
							}
						}
					}
				}
				else {
					$(built.node).append($('<span class="dt-down-arrow">'+this.c.dom.splitDropdown.text+'</span>'))
				}

				this._expandButton( built.buttons, built.conf.buttons, built.conf.split, !isSplit, isSplit, attachPoint, built.conf );
			}
			built.conf.parent = parentConf;

			// init call is made here, rather than buildButton as it needs to
			// be selectable, and for that it needs to be in the buttons array
			if ( conf.init ) {
				conf.init.call( dt.button( built.node ), dt, $(built.node), conf );
			}

			buttonCounter++;
		}
	},

	/**
	 * Create an individual button
	 * @param  {object} config            Resolved button configuration
	 * @param  {boolean} inCollection `true` if a collection button
	 * @return {jQuery} Created button node (jQuery)
	 * @private
	 */
	_buildButton: function ( config, inCollection, isSplit, inSplit )
	{
		var buttonDom = this.c.dom.button;
		var linerDom = this.c.dom.buttonLiner;
		var collectionDom = this.c.dom.collection;
		var splitDom = this.c.dom.split;
		var splitCollectionDom = this.c.dom.splitCollection;
		var splitDropdownButton = this.c.dom.splitDropdownButton;
		var dt = this.s.dt;
		var text = function ( opt ) {
			return typeof opt === 'function' ?
				opt( dt, button, config ) :
				opt;
		};

		// Spacers don't do much other than insert an element into the DOM
		if (config.spacer) {
			var spacer = $('<span></span>')
				.addClass('dt-button-spacer ' + config.style + ' ' + buttonDom.spacerClass)
				.html(text(config.text));

			return {
				conf:         config,
				node:         spacer,
				inserter:     spacer,
				buttons:      [],
				inCollection: inCollection,
				isSplit:	  isSplit,
				inSplit:	  inSplit,
				collection:   null
			};
		}

		if ( !isSplit && inSplit && splitCollectionDom ) {
			buttonDom = splitDropdownButton;
		}
		else if ( !isSplit && inCollection && collectionDom.button ) {
			buttonDom = collectionDom.button;
		} 

		if ( !isSplit && inSplit && splitCollectionDom.buttonLiner ) {
			linerDom = splitCollectionDom.buttonLiner
		}
		else if ( !isSplit && inCollection && collectionDom.buttonLiner ) {
			linerDom = collectionDom.buttonLiner;
		}

		// Make sure that the button is available based on whatever requirements
		// it has. For example, PDF button require pdfmake
		if ( config.available && ! config.available( dt, config ) && !config.hasOwnProperty('html') ) {
			return false;
		}

		var button;
		if(!config.hasOwnProperty('html')) {
			var action = function ( e, dt, button, config ) {
				config.action.call( dt.button( button ), e, dt, button, config );
	
				$(dt.table().node()).triggerHandler( 'buttons-action.dt', [
					dt.button( button ), dt, button, config 
				] );
			};

			var tag = config.tag || buttonDom.tag;
			var clickBlurs = config.clickBlurs === undefined
				? true :
				config.clickBlurs;

			button = $('<'+tag+'/>')
				.addClass( buttonDom.className )
				.addClass( inSplit ? this.c.dom.splitDropdownButton.className : '')
				.attr( 'tabindex', this.s.dt.settings()[0].iTabIndex )
				.attr( 'aria-controls', this.s.dt.table().node().id )
				.on( 'click.dtb', function (e) {
					e.preventDefault();
	
					if ( ! button.hasClass( buttonDom.disabled ) && config.action ) {
						action( e, dt, button, config );
					}
					if( clickBlurs ) {
						button.trigger('blur');
					}
				} )
				.on( 'keypress.dtb', function (e) {
					if ( e.keyCode === 13 ) {
						e.preventDefault();

						if ( ! button.hasClass( buttonDom.disabled ) && config.action ) {
							action( e, dt, button, config );
						}
					}
				} );
	
			// Make `a` tags act like a link
			if ( tag.toLowerCase() === 'a' ) {
				button.attr( 'href', '#' );
			}
	
			// Button tags should have `type=button` so they don't have any default behaviour
			if ( tag.toLowerCase() === 'button' ) {
				button.attr( 'type', 'button' );
			}
	
			if ( linerDom.tag ) {
				var liner = $('<'+linerDom.tag+'/>')
					.html( text( config.text ) )
					.addClass( linerDom.className );
	
				if ( linerDom.tag.toLowerCase() === 'a' ) {
					liner.attr( 'href', '#' );
				}
	
				button.append( liner );
			}
			else {
				button.html( text( config.text ) );
			}
	
			if ( config.enabled === false ) {
				button.addClass( buttonDom.disabled );
			}
	
			if ( config.className ) {
				button.addClass( config.className );
			}
	
			if ( config.titleAttr ) {
				button.attr( 'title', text( config.titleAttr ) );
			}
	
			if ( config.attr ) {
				button.attr( config.attr );
			}
	
			if ( ! config.namespace ) {
				config.namespace = '.dt-button-'+(_buttonCounter++);
			}

			if  ( config.config !== undefined && config.config.split ) {
				config.split = config.config.split;
			}
		}
		else {
			button = $(config.html)
		}
	
		var buttonContainer = this.c.dom.buttonContainer;
		var inserter;
		if ( buttonContainer && buttonContainer.tag ) {
			inserter = $('<'+buttonContainer.tag+'/>')
				.addClass( buttonContainer.className )
				.append( button );
		}
		else {
			inserter = button;
		}

		this._addKey( config );

		// Style integration callback for DOM manipulation
		// Note that this is _not_ documented. It is currently
		// for style integration only
		if( this.c.buttonCreated ) {
			inserter = this.c.buttonCreated( config, inserter );
		}

		var splitDiv;
		if(isSplit) {
			splitDiv = $('<div/>').addClass(this.c.dom.splitWrapper.className)
			splitDiv.append(button);
			var dropButtonConfig = $.extend(config, {
				text: this.c.dom.splitDropdown.text,
				className: this.c.dom.splitDropdown.className,
				closeButton: false,
				attr: {
					'aria-haspopup': 'dialog',
					'aria-expanded': false
				},
				align: this.c.dom.splitDropdown.align,
				splitAlignClass: this.c.dom.splitDropdown.splitAlignClass
				
			})

			this._addKey(dropButtonConfig);

			var splitAction = function ( e, dt, button, config ) {
				_dtButtons.split.action.call( dt.button(splitDiv), e, dt, button, config );
	
				$(dt.table().node()).triggerHandler( 'buttons-action.dt', [
					dt.button( button ), dt, button, config 
				] );
				button.attr('aria-expanded', true)
			};
			
			var dropButton = $('<button class="' + this.c.dom.splitDropdown.className + ' dt-button"><span class="dt-btn-split-drop-arrow">'+this.c.dom.splitDropdown.text+'</span></button>')
				.on( 'click.dtb', function (e) {
					e.preventDefault();
					e.stopPropagation();

					if ( ! dropButton.hasClass( buttonDom.disabled )) {
						splitAction( e, dt, dropButton, dropButtonConfig );
					}
					if ( clickBlurs ) {
						dropButton.trigger('blur');
					}
				} )
				.on( 'keypress.dtb', function (e) {
					if ( e.keyCode === 13 ) {
						e.preventDefault();

						if ( ! dropButton.hasClass( buttonDom.disabled ) ) {
							splitAction( e, dt, dropButton, dropButtonConfig );
						}
					}
				} );

			if(config.split.length === 0) {
				dropButton.addClass('dtb-hide-drop');
			}

			splitDiv.append(dropButton).attr(dropButtonConfig.attr);
		}

		return {
			conf:         config,
			node:         isSplit ? splitDiv.get(0) : button.get(0),
			inserter:     isSplit ? splitDiv : inserter,
			buttons:      [],
			inCollection: inCollection,
			isSplit:	  isSplit,
			inSplit:	  inSplit,
			collection:   null
		};
	},

	/**
	 * Get the button object from a node (recursive)
	 * @param  {node} node Button node
	 * @param  {array} [buttons] Button array, uses base if not defined
	 * @return {object} Button object
	 * @private
	 */
	_nodeToButton: function ( node, buttons )
	{
		if ( ! buttons ) {
			buttons = this.s.buttons;
		}

		for ( var i=0, ien=buttons.length ; i<ien ; i++ ) {
			if ( buttons[i].node === node ) {
				return buttons[i];
			}

			if ( buttons[i].buttons.length ) {
				var ret = this._nodeToButton( node, buttons[i].buttons );

				if ( ret ) {
					return ret;
				}
			}
		}
	},

	/**
	 * Get container array for a button from a button node (recursive)
	 * @param  {node} node Button node
	 * @param  {array} [buttons] Button array, uses base if not defined
	 * @return {array} Button's host array
	 * @private
	 */
	_nodeToHost: function ( node, buttons )
	{
		if ( ! buttons ) {
			buttons = this.s.buttons;
		}

		for ( var i=0, ien=buttons.length ; i<ien ; i++ ) {
			if ( buttons[i].node === node ) {
				return buttons;
			}

			if ( buttons[i].buttons.length ) {
				var ret = this._nodeToHost( node, buttons[i].buttons );

				if ( ret ) {
					return ret;
				}
			}
		}
	},

	/**
	 * Handle a key press - determine if any button's key configured matches
	 * what was typed and trigger the action if so.
	 * @param  {string} character The character pressed
	 * @param  {object} e Key event that triggered this call
	 * @private
	 */
	_keypress: function ( character, e )
	{
		// Check if this button press already activated on another instance of Buttons
		if ( e._buttonsHandled ) {
			return;
		}

		var run = function ( conf, node ) {
			if ( ! conf.key ) {
				return;
			}

			if ( conf.key === character ) {
				e._buttonsHandled = true;
				$(node).click();
			}
			else if ( $.isPlainObject( conf.key ) ) {
				if ( conf.key.key !== character ) {
					return;
				}

				if ( conf.key.shiftKey && ! e.shiftKey ) {
					return;
				}

				if ( conf.key.altKey && ! e.altKey ) {
					return;
				}

				if ( conf.key.ctrlKey && ! e.ctrlKey ) {
					return;
				}

				if ( conf.key.metaKey && ! e.metaKey ) {
					return;
				}

				// Made it this far - it is good
				e._buttonsHandled = true;
				$(node).click();
			}
		};

		var recurse = function ( a ) {
			for ( var i=0, ien=a.length ; i<ien ; i++ ) {
				run( a[i].conf, a[i].node );

				if ( a[i].buttons.length ) {
					recurse( a[i].buttons );
				}
			}
		};

		recurse( this.s.buttons );
	},

	/**
	 * Remove a key from the key listener for this instance (to be used when a
	 * button is removed)
	 * @param  {object} conf Button configuration
	 * @private
	 */
	_removeKey: function ( conf )
	{
		if ( conf.key ) {
			var character = $.isPlainObject( conf.key ) ?
				conf.key.key :
				conf.key;

			// Remove only one character, as multiple buttons could have the
			// same listening key
			var a = this.s.listenKeys.split('');
			var idx = $.inArray( character, a );
			a.splice( idx, 1 );
			this.s.listenKeys = a.join('');
		}
	},

	/**
	 * Resolve a button configuration
	 * @param  {string|function|object} conf Button config to resolve
	 * @return {object} Button configuration
	 * @private
	 */
	_resolveExtends: function ( conf )
	{
		var that = this;
		var dt = this.s.dt;
		var i, ien;
		var toConfObject = function ( base ) {
			var loop = 0;

			// Loop until we have resolved to a button configuration, or an
			// array of button configurations (which will be iterated
			// separately)
			while ( ! $.isPlainObject(base) && ! Array.isArray(base) ) {
				if ( base === undefined ) {
					return;
				}

				if ( typeof base === 'function' ) {
					base = base.call( that, dt, conf );

					if ( ! base ) {
						return false;
					}
				}
				else if ( typeof base === 'string' ) {
					if ( ! _dtButtons[ base ] ) {
						return {html: base}
					}

					base = _dtButtons[ base ];
				}

				loop++;
				if ( loop > 30 ) {
					// Protect against misconfiguration killing the browser
					throw 'Buttons: Too many iterations';
				}
			}

			return Array.isArray( base ) ?
				base :
				$.extend( {}, base );
		};

		conf = toConfObject( conf );

		while ( conf && conf.extend ) {
			// Use `toConfObject` in case the button definition being extended
			// is itself a string or a function
			if ( ! _dtButtons[ conf.extend ] ) {
				throw 'Cannot extend unknown button type: '+conf.extend;
			}

			var objArray = toConfObject( _dtButtons[ conf.extend ] );
			if ( Array.isArray( objArray ) ) {
				return objArray;
			}
			else if ( ! objArray ) {
				// This is a little brutal as it might be possible to have a
				// valid button without the extend, but if there is no extend
				// then the host button would be acting in an undefined state
				return false;
			}

			// Stash the current class name
			var originalClassName = objArray.className;

			if (conf.config !== undefined && objArray.config !== undefined) {
				conf.config = $.extend({}, objArray.config, conf.config)
			}

			conf = $.extend( {}, objArray, conf );

			// The extend will have overwritten the original class name if the
			// `conf` object also assigned a class, but we want to concatenate
			// them so they are list that is combined from all extended buttons
			if ( originalClassName && conf.className !== originalClassName ) {
				conf.className = originalClassName+' '+conf.className;
			}

			// Although we want the `conf` object to overwrite almost all of
			// the properties of the object being extended, the `extend`
			// property should come from the object being extended
			conf.extend = objArray.extend;
		}

		// Buttons to be added to a collection  -gives the ability to define
		// if buttons should be added to the start or end of a collection
		var postfixButtons = conf.postfixButtons;
		if ( postfixButtons ) {
			if ( ! conf.buttons ) {
				conf.buttons = [];
			}

			for ( i=0, ien=postfixButtons.length ; i<ien ; i++ ) {
				conf.buttons.push( postfixButtons[i] );
			}
		}

		var prefixButtons = conf.prefixButtons;
		if ( prefixButtons ) {
			if ( ! conf.buttons ) {
				conf.buttons = [];
			}

			for ( i=0, ien=prefixButtons.length ; i<ien ; i++ ) {
				conf.buttons.splice( i, 0, prefixButtons[i] );
			}
		}

		return conf;
	},

	/**
	 * Display (and replace if there is an existing one) a popover attached to a button
	 * @param {string|node} content Content to show
	 * @param {DataTable.Api} hostButton DT API instance of the button
	 * @param {object} inOpts Options (see object below for all options)
	 */
	_popover: function ( content, hostButton, inOpts, e ) {
		var dt = hostButton;
		var buttonsSettings = this.c;
		var closed = false;
		var options = $.extend( {
			align: 'button-left', // button-right, dt-container, split-left, split-right
			autoClose: false,
			background: true,
			backgroundClassName: 'dt-button-background',
			closeButton: true,
			contentClassName: buttonsSettings.dom.collection.className,
			collectionLayout: '',
			collectionTitle: '',
			dropup: false,
			fade: 400,
			popoverTitle: '',
			rightAlignClassName: 'dt-button-right',
			tag: buttonsSettings.dom.collection.tag
		}, inOpts );

		var hostNode = hostButton.node();

		var close = function () {
			closed = true;

			_fadeOut(
				$('.dt-button-collection'),
				options.fade,
				function () {
					$(this).detach();
				}
			);

			$(dt.buttons( '[aria-haspopup="dialog"][aria-expanded="true"]' ).nodes())
				.attr('aria-expanded', 'false');

			$('div.dt-button-background').off( 'click.dtb-collection' );
			Buttons.background( false, options.backgroundClassName, options.fade, hostNode );

			$(window).off('resize.resize.dtb-collection');
			$('body').off( '.dtb-collection' );
			dt.off( 'buttons-action.b-internal' );
			dt.off( 'destroy' );
		};

		if (content === false) {
			close();
			return;
		}

		var existingExpanded = $(dt.buttons( '[aria-haspopup="dialog"][aria-expanded="true"]' ).nodes());
		if ( existingExpanded.length ) {
			// Reuse the current position if the button that was triggered is inside an existing collection
			if (hostNode.closest('div.dt-button-collection').length) {
				hostNode = existingExpanded.eq(0);
			}

			close();
		}

		// Try to be smart about the layout
		var cnt = $('.dt-button', content).length;
		var mod = '';

		if (cnt === 3) {
			mod = 'dtb-b3';
		}
		else if (cnt === 2) {
			mod = 'dtb-b2';
		}
		else if (cnt === 1) {
			mod = 'dtb-b1';
		}

		var display = $('<div/>')
			.addClass('dt-button-collection')
			.addClass(options.collectionLayout)
			.addClass(options.splitAlignClass)
			.addClass(mod)
			.css('display', 'none')
			.attr({
				'aria-modal': true,
				role: 'dialog'
			});

		content = $(content)
			.addClass(options.contentClassName)
			.attr('role', 'menu')
			.appendTo(display);

		hostNode.attr( 'aria-expanded', 'true' );

		if ( hostNode.parents('body')[0] !== document.body ) {
			hostNode = document.body.lastChild;
		}

		if ( options.popoverTitle ) {
			display.prepend('<div class="dt-button-collection-title">'+options.popoverTitle+'</div>');
		}
		else if ( options.collectionTitle ) {
			display.prepend('<div class="dt-button-collection-title">'+options.collectionTitle+'</div>');
		}

		if (options.closeButton) {
			display.prepend('<div class="dtb-popover-close">x</div>').addClass('dtb-collection-closeable')
		}

		_fadeIn( display.insertAfter( hostNode ), options.fade );

		var tableContainer = $( hostButton.table().container() );
		var position = display.css( 'position' );

		if ( options.span === 'container' || options.align === 'dt-container' ) {
			hostNode = hostNode.parent();
			display.css('width', tableContainer.width());
		}

		// Align the popover relative to the DataTables container
		// Useful for wide popovers such as SearchPanes
		if (position === 'absolute') {
			// Align relative to the host button
			var offsetParent = $(hostNode[0].offsetParent);
			var buttonPosition = hostNode.position();
			var buttonOffset = hostNode.offset();
			var tableSizes = offsetParent.offset();
			var containerPosition = offsetParent.position();
			var computed = window.getComputedStyle(offsetParent[0]);

			tableSizes.height = offsetParent.outerHeight();
			tableSizes.width = offsetParent.width() + parseFloat(computed.paddingLeft);
			tableSizes.right = tableSizes.left + tableSizes.width;
			tableSizes.bottom = tableSizes.top + tableSizes.height;

			// Set the initial position so we can read height / width
			var top = buttonPosition.top + hostNode.outerHeight();
			var left = buttonPosition.left;

			display.css( {
				top: top,
				left: left
			} );

			// Get the popover position
			computed = window.getComputedStyle(display[0]);
			var popoverSizes = display.offset();

			popoverSizes.height = display.outerHeight();
			popoverSizes.width = display.outerWidth();
			popoverSizes.right = popoverSizes.left + popoverSizes.width;
			popoverSizes.bottom = popoverSizes.top + popoverSizes.height;
			popoverSizes.marginTop = parseFloat(computed.marginTop);
			popoverSizes.marginBottom = parseFloat(computed.marginBottom);

			// First position per the class requirements - pop up and right align
			if (options.dropup) {
				top = buttonPosition.top - popoverSizes.height - popoverSizes.marginTop - popoverSizes.marginBottom;
			}

			if (options.align === 'button-right' || display.hasClass( options.rightAlignClassName )) {
				left = buttonPosition.left - popoverSizes.width + hostNode.outerWidth(); 
			}

			// Container alignment - make sure it doesn't overflow the table container
			if (options.align === 'dt-container' || options.align === 'container') {
				if (left < buttonPosition.left) {
					left = -buttonPosition.left;
				}

				if (left + popoverSizes.width > tableSizes.width) {
					left = tableSizes.width - popoverSizes.width;
				}
			}

			// Window adjustment
			if (containerPosition.left + left + popoverSizes.width > $(window).width()) {
				// Overflowing the document to the right
				left = $(window).width() - popoverSizes.width - containerPosition.left;
			}

			if (buttonOffset.left + left < 0) {
				// Off to the left of the document
				left = -buttonOffset.left;
			}

			if (containerPosition.top + top + popoverSizes.height > $(window).height() + $(window).scrollTop()) {
				// Pop up if otherwise we'd need the user to scroll down
				top = buttonPosition.top - popoverSizes.height - popoverSizes.marginTop - popoverSizes.marginBottom;
			}

			if (containerPosition.top + top < $(window).scrollTop()) {
				// Correction for when the top is beyond the top of the page
				top = buttonPosition.top + hostNode.outerHeight();
			}

			// Calculations all done - now set it
			display.css( {
				top: top,
				left: left
			} );
		}
		else {
			// Fix position - centre on screen
			var position = function () {
				var half = $(window).height() / 2;

				var top = display.height() / 2;
				if ( top > half ) {
					top = half;
				}

				display.css( 'marginTop', top*-1 );
			};

			position();

			$(window).on('resize.dtb-collection', function () {
				position();
			});
		}

		if ( options.background ) {
			Buttons.background(
				true,
				options.backgroundClassName,
				options.fade,
				options.backgroundHost || hostNode
			);
		}

		// This is bonkers, but if we don't have a click listener on the
		// background element, iOS Safari will ignore the body click
		// listener below. An empty function here is all that is
		// required to make it work...
		$('div.dt-button-background').on( 'click.dtb-collection', function () {} );

		if ( options.autoClose ) {
			setTimeout( function () {
				dt.on( 'buttons-action.b-internal', function (e, btn, dt, node) {
					if ( node[0] === hostNode[0] ) {
						return;
					}
					close();
				} );
			}, 0);
		}
		
		$(display).trigger('buttons-popover.dt');


		dt.on('destroy', close);

		setTimeout(function() {
			closed = false;
			$('body')
				.on( 'click.dtb-collection', function (e) {
					if (closed) {
						return;
					}

					// andSelf is deprecated in jQ1.8, but we want 1.7 compat
					var back = $.fn.addBack ? 'addBack' : 'andSelf';
					var parent = $(e.target).parent()[0];
	
					if (( ! $(e.target).parents()[back]().filter( content ).length  && !$(parent).hasClass('dt-buttons')) || $(e.target).hasClass('dt-button-background')) {
						close();
					}
				} )
				.on( 'keyup.dtb-collection', function (e) {
					if ( e.keyCode === 27 ) {
						close();
					}
				} )
				.on( 'keydown.dtb-collection', function (e) {
					// Focus trap for tab key
					var elements = $('a, button', content);
					var active = document.activeElement;

					if (e.keyCode !== 9) { // tab
						return;
					}

					if (elements.index(active) === -1) {
						// If current focus is not inside the popover
						elements.first().focus();
						e.preventDefault();
					}
					else if (e.shiftKey) {
						// Reverse tabbing order when shift key is pressed
						if (active === elements[0]) {
							elements.last().focus();
							e.preventDefault();
						}
					}
					else {
						if (active === elements.last()[0]) {
							elements.first().focus();
							e.preventDefault();
						}
					}
				} );
		}, 0);
	}
} );



/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Statics
 */

/**
 * Show / hide a background layer behind a collection
 * @param  {boolean} Flag to indicate if the background should be shown or
 *   hidden 
 * @param  {string} Class to assign to the background
 * @static
 */
Buttons.background = function ( show, className, fade, insertPoint ) {
	if ( fade === undefined ) {
		fade = 400;
	}
	if ( ! insertPoint ) {
		insertPoint = document.body;
	}

	if ( show ) {
		_fadeIn(
			$('<div/>')
				.addClass( className )
				.css( 'display', 'none' )
				.insertAfter( insertPoint ),
			fade
		);
	}
	else {
		_fadeOut(
			$('div.'+className),
			fade,
			function () {
				$(this)
					.removeClass( className )
					.remove();
			}
		);
	}
};

/**
 * Instance selector - select Buttons instances based on an instance selector
 * value from the buttons assigned to a DataTable. This is only useful if
 * multiple instances are attached to a DataTable.
 * @param  {string|int|array} Instance selector - see `instance-selector`
 *   documentation on the DataTables site
 * @param  {array} Button instance array that was attached to the DataTables
 *   settings object
 * @return {array} Buttons instances
 * @static
 */
Buttons.instanceSelector = function ( group, buttons )
{
	if ( group === undefined || group === null ) {
		return $.map( buttons, function ( v ) {
			return v.inst;
		} );
	}

	var ret = [];
	var names = $.map( buttons, function ( v ) {
		return v.name;
	} );

	// Flatten the group selector into an array of single options
	var process = function ( input ) {
		if ( Array.isArray( input ) ) {
			for ( var i=0, ien=input.length ; i<ien ; i++ ) {
				process( input[i] );
			}
			return;
		}

		if ( typeof input === 'string' ) {
			if ( input.indexOf( ',' ) !== -1 ) {
				// String selector, list of names
				process( input.split(',') );
			}
			else {
				// String selector individual name
				var idx = $.inArray( input.trim(), names );

				if ( idx !== -1 ) {
					ret.push( buttons[ idx ].inst );
				}
			}
		}
		else if ( typeof input === 'number' ) {
			// Index selector
			ret.push( buttons[ input ].inst );
		}
		else if ( typeof input === 'object' ) {
			// Actual instance selector
			ret.push( input );
		}
	};
	
	process( group );

	return ret;
};

/**
 * Button selector - select one or more buttons from a selector input so some
 * operation can be performed on them.
 * @param  {array} Button instances array that the selector should operate on
 * @param  {string|int|node|jQuery|array} Button selector - see
 *   `button-selector` documentation on the DataTables site
 * @return {array} Array of objects containing `inst` and `idx` properties of
 *   the selected buttons so you know which instance each button belongs to.
 * @static
 */
Buttons.buttonSelector = function ( insts, selector )
{
	var ret = [];
	var nodeBuilder = function ( a, buttons, baseIdx ) {
		var button;
		var idx;

		for ( var i=0, ien=buttons.length ; i<ien ; i++ ) {
			button = buttons[i];

			if ( button ) {
				idx = baseIdx !== undefined ?
					baseIdx+i :
					i+'';

				a.push( {
					node: button.node,
					name: button.conf.name,
					idx:  idx
				} );

				if ( button.buttons ) {
					nodeBuilder( a, button.buttons, idx+'-' );
				}
			}
		}
	};

	var run = function ( selector, inst ) {
		var i, ien;
		var buttons = [];
		nodeBuilder( buttons, inst.s.buttons );

		var nodes = $.map( buttons, function (v) {
			return v.node;
		} );

		if ( Array.isArray( selector ) || selector instanceof $ ) {
			for ( i=0, ien=selector.length ; i<ien ; i++ ) {
				run( selector[i], inst );
			}
			return;
		}

		if ( selector === null || selector === undefined || selector === '*' ) {
			// Select all
			for ( i=0, ien=buttons.length ; i<ien ; i++ ) {
				ret.push( {
					inst: inst,
					node: buttons[i].node
				} );
			}
		}
		else if ( typeof selector === 'number' ) {
			// Main button index selector
			if (inst.s.buttons[ selector ]) {
				ret.push( {
					inst: inst,
					node: inst.s.buttons[ selector ].node
				} );
			}
		}
		else if ( typeof selector === 'string' ) {
			if ( selector.indexOf( ',' ) !== -1 ) {
				// Split
				var a = selector.split(',');

				for ( i=0, ien=a.length ; i<ien ; i++ ) {
					run( a[i].trim(), inst );
				}
			}
			else if ( selector.match( /^\d+(\-\d+)*$/ ) ) {
				// Sub-button index selector
				var indexes = $.map( buttons, function (v) {
					return v.idx;
				} );

				ret.push( {
					inst: inst,
					node: buttons[ $.inArray( selector, indexes ) ].node
				} );
			}
			else if ( selector.indexOf( ':name' ) !== -1 ) {
				// Button name selector
				var name = selector.replace( ':name', '' );

				for ( i=0, ien=buttons.length ; i<ien ; i++ ) {
					if ( buttons[i].name === name ) {
						ret.push( {
							inst: inst,
							node: buttons[i].node
						} );
					}
				}
			}
			else {
				// jQuery selector on the nodes
				$( nodes ).filter( selector ).each( function () {
					ret.push( {
						inst: inst,
						node: this
					} );
				} );
			}
		}
		else if ( typeof selector === 'object' && selector.nodeName ) {
			// Node selector
			var idx = $.inArray( selector, nodes );

			if ( idx !== -1 ) {
				ret.push( {
					inst: inst,
					node: nodes[ idx ]
				} );
			}
		}
	};


	for ( var i=0, ien=insts.length ; i<ien ; i++ ) {
		var inst = insts[i];

		run( selector, inst );
	}

	return ret;
};

/**
 * Default function used for formatting output data.
 * @param {*} str Data to strip
 */
Buttons.stripData = function ( str, config ) {
	if ( typeof str !== 'string' ) {
		return str;
	}

	// Always remove script tags
	str = str.replace( /<script\b[^<]*(?:(?!<\/script>)<[^<]*)*<\/script>/gi, '' );

	// Always remove comments
	str = str.replace( /<!\-\-.*?\-\->/g, '' );

	if ( ! config || config.stripHtml ) {
		str = str.replace( /<[^>]*>/g, '' );
	}

	if ( ! config || config.trim ) {
		str = str.replace( /^\s+|\s+$/g, '' );
	}

	if ( ! config || config.stripNewlines ) {
		str = str.replace( /\n/g, ' ' );
	}

	if ( ! config || config.decodeEntities ) {
		_exportTextarea.innerHTML = str;
		str = _exportTextarea.value;
	}

	return str;
};


/**
 * Buttons defaults. For full documentation, please refer to the docs/option
 * directory or the DataTables site.
 * @type {Object}
 * @static
 */
Buttons.defaults = {
	buttons: [ 'copy', 'excel', 'csv', 'pdf', 'print' ],
	name: 'main',
	tabIndex: 0,
	dom: {
		container: {
			tag: 'div',
			className: 'dt-buttons'
		},
		collection: {
			tag: 'div',
			className: ''
		},
		button: {
			tag: 'button',
			className: 'dt-button',
			active: 'active',
			disabled: 'disabled',
			spacerClass: ''
		},
		buttonLiner: {
			tag: 'span',
			className: ''
		},
		split: {
			tag: 'div',
			className: 'dt-button-split',
		},
		splitWrapper: {
			tag: 'div',
			className: 'dt-btn-split-wrapper',
		},
		splitDropdown: {
			tag: 'button',
			text: '&#x25BC;',
			className: 'dt-btn-split-drop',
			align: 'split-right',
			splitAlignClass: 'dt-button-split-left'
		},
		splitDropdownButton: {
			tag: 'button',
			className: 'dt-btn-split-drop-button dt-button',
		},
		splitCollection: {
			tag: 'div',
			className: 'dt-button-split-collection',
		}
	}
};

/**
 * Version information
 * @type {string}
 * @static
 */
Buttons.version = '2.3.6';


$.extend( _dtButtons, {
	collection: {
		text: function ( dt ) {
			return dt.i18n( 'buttons.collection', 'Collection' );
		},
		className: 'buttons-collection',
		closeButton: false,
		init: function ( dt, button, config ) {
			button.attr( 'aria-expanded', false );
		},
		action: function ( e, dt, button, config ) {
			if ( config._collection.parents('body').length ) {
				this.popover(false, config);
			}
			else {
				this.popover(config._collection, config);
			}

			// When activated using a key - auto focus on the
			// first item in the popover
			if (e.type === 'keypress') {
				$('a, button', config._collection).eq(0).focus();
			}
		},
		attr: {
			'aria-haspopup': 'dialog'
		}
		// Also the popover options, defined in Buttons.popover
	},
	split: {
		text: function ( dt ) {
			return dt.i18n( 'buttons.split', 'Split' );
		},
		className: 'buttons-split',
		closeButton: false,
		init: function ( dt, button, config ) {
			return button.attr( 'aria-expanded', false );
		},
		action: function ( e, dt, button, config ) {
			this.popover(config._collection, config);
		},
		attr: {
			'aria-haspopup': 'dialog'
		}
		// Also the popover options, defined in Buttons.popover
	},
	copy: function ( dt, conf ) {
		if ( _dtButtons.copyHtml5 ) {
			return 'copyHtml5';
		}
	},
	csv: function ( dt, conf ) {
		if ( _dtButtons.csvHtml5 && _dtButtons.csvHtml5.available( dt, conf ) ) {
			return 'csvHtml5';
		}
	},
	excel: function ( dt, conf ) {
		if ( _dtButtons.excelHtml5 && _dtButtons.excelHtml5.available( dt, conf ) ) {
			return 'excelHtml5';
		}
	},
	pdf: function ( dt, conf ) {
		if ( _dtButtons.pdfHtml5 && _dtButtons.pdfHtml5.available( dt, conf ) ) {
			return 'pdfHtml5';
		}
	},
	pageLength: function ( dt ) {
		var lengthMenu = dt.settings()[0].aLengthMenu;
		var vals = [];
		var lang = [];
		var text = function ( dt ) {
			return dt.i18n( 'buttons.pageLength', {
				"-1": 'Show all rows',
				_:    'Show %d rows'
			}, dt.page.len() );
		};

		// Support for DataTables 1.x 2D array
		if (Array.isArray( lengthMenu[0] )) {
			vals = lengthMenu[0];
			lang = lengthMenu[1];
		}
		else {
			for (var i=0 ; i<lengthMenu.length ; i++) {
				var option = lengthMenu[i];

				// Support for DataTables 2 object in the array
				if ($.isPlainObject(option)) {
					vals.push(option.value);
					lang.push(option.label);
				}
				else {
					vals.push(option);
					lang.push(option);
				}
			}
		}

		return {
			extend: 'collection',
			text: text,
			className: 'buttons-page-length',
			autoClose: true,
			buttons: $.map( vals, function ( val, i ) {
				return {
					text: lang[i],
					className: 'button-page-length',
					action: function ( e, dt ) {
						dt.page.len( val ).draw();
					},
					init: function ( dt, node, conf ) {
						var that = this;
						var fn = function () {
							that.active( dt.page.len() === val );
						};

						dt.on( 'length.dt'+conf.namespace, fn );
						fn();
					},
					destroy: function ( dt, node, conf ) {
						dt.off( 'length.dt'+conf.namespace );
					}
				};
			} ),
			init: function ( dt, node, conf ) {
				var that = this;
				dt.on( 'length.dt'+conf.namespace, function () {
					that.text( conf.text );
				} );
			},
			destroy: function ( dt, node, conf ) {
				dt.off( 'length.dt'+conf.namespace );
			}
		};
	},
	spacer: {
		style: 'empty',
		spacer: true,
		text: function ( dt ) {
			return dt.i18n( 'buttons.spacer', '' );
		}
	}
} );


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * DataTables API
 *
 * For complete documentation, please refer to the docs/api directory or the
 * DataTables site
 */

// Buttons group and individual button selector
DataTable.Api.register( 'buttons()', function ( group, selector ) {
	// Argument shifting
	if ( selector === undefined ) {
		selector = group;
		group = undefined;
	}

	this.selector.buttonGroup = group;

	var res = this.iterator( true, 'table', function ( ctx ) {
		if ( ctx._buttons ) {
			return Buttons.buttonSelector(
				Buttons.instanceSelector( group, ctx._buttons ),
				selector
			);
		}
	}, true );

	res._groupSelector = group;
	return res;
} );

// Individual button selector
DataTable.Api.register( 'button()', function ( group, selector ) {
	// just run buttons() and truncate
	var buttons = this.buttons( group, selector );

	if ( buttons.length > 1 ) {
		buttons.splice( 1, buttons.length );
	}

	return buttons;
} );

// Active buttons
DataTable.Api.registerPlural( 'buttons().active()', 'button().active()', function ( flag ) {
	if ( flag === undefined ) {
		return this.map( function ( set ) {
			return set.inst.active( set.node );
		} );
	}

	return this.each( function ( set ) {
		set.inst.active( set.node, flag );
	} );
} );

// Get / set button action
DataTable.Api.registerPlural( 'buttons().action()', 'button().action()', function ( action ) {
	if ( action === undefined ) {
		return this.map( function ( set ) {
			return set.inst.action( set.node );
		} );
	}

	return this.each( function ( set ) {
		set.inst.action( set.node, action );
	} );
} );

// Collection control
DataTable.Api.registerPlural( 'buttons().collectionRebuild()', 'button().collectionRebuild()', function ( buttons ) {
	return this.each( function ( set ) {
		for(var i = 0; i < buttons.length; i++) {
			if(typeof buttons[i] === 'object') {
				buttons[i].parentConf = set;
			}
		}
		set.inst.collectionRebuild( set.node, buttons );
	} );
} );

// Enable / disable buttons
DataTable.Api.register( ['buttons().enable()', 'button().enable()'], function ( flag ) {
	return this.each( function ( set ) {
		set.inst.enable( set.node, flag );
	} );
} );

// Disable buttons
DataTable.Api.register( ['buttons().disable()', 'button().disable()'], function () {
	return this.each( function ( set ) {
		set.inst.disable( set.node );
	} );
} );

// Button index
DataTable.Api.register( 'button().index()', function () {
	var idx = null;

	this.each( function ( set ) {
		var res = set.inst.index( set.node );

		if (res !== null) {
			idx = res;
		}
	} );

	return idx;
} );

// Get button nodes
DataTable.Api.registerPlural( 'buttons().nodes()', 'button().node()', function () {
	var jq = $();

	// jQuery will automatically reduce duplicates to a single entry
	$( this.each( function ( set ) {
		jq = jq.add( set.inst.node( set.node ) );
	} ) );

	return jq;
} );

// Get / set button processing state
DataTable.Api.registerPlural( 'buttons().processing()', 'button().processing()', function ( flag ) {
	if ( flag === undefined ) {
		return this.map( function ( set ) {
			return set.inst.processing( set.node );
		} );
	}

	return this.each( function ( set ) {
		set.inst.processing( set.node, flag );
	} );
} );

// Get / set button text (i.e. the button labels)
DataTable.Api.registerPlural( 'buttons().text()', 'button().text()', function ( label ) {
	if ( label === undefined ) {
		return this.map( function ( set ) {
			return set.inst.text( set.node );
		} );
	}

	return this.each( function ( set ) {
		set.inst.text( set.node, label );
	} );
} );

// Trigger a button's action
DataTable.Api.registerPlural( 'buttons().trigger()', 'button().trigger()', function () {
	return this.each( function ( set ) {
		set.inst.node( set.node ).trigger( 'click' );
	} );
} );

// Button resolver to the popover
DataTable.Api.register( 'button().popover()', function (content, options) {
	return this.map( function ( set ) {
		return set.inst._popover( content, this.button(this[0].node), options );
	} );
} );

// Get the container elements
DataTable.Api.register( 'buttons().containers()', function () {
	var jq = $();
	var groupSelector = this._groupSelector;

	// We need to use the group selector directly, since if there are no buttons
	// the result set will be empty
	this.iterator( true, 'table', function ( ctx ) {
		if ( ctx._buttons ) {
			var insts = Buttons.instanceSelector( groupSelector, ctx._buttons );

			for ( var i=0, ien=insts.length ; i<ien ; i++ ) {
				jq = jq.add( insts[i].container() );
			}
		}
	} );

	return jq;
} );

DataTable.Api.register( 'buttons().container()', function () {
	// API level of nesting is `buttons()` so we can zip into the containers method
	return this.containers().eq(0);
} );

// Add a new button
DataTable.Api.register( 'button().add()', function ( idx, conf, draw ) {
	var ctx = this.context;

	// Don't use `this` as it could be empty - select the instances directly
	if ( ctx.length ) {
		var inst = Buttons.instanceSelector( this._groupSelector, ctx[0]._buttons );

		if ( inst.length ) {
			inst[0].add( conf, idx , draw);
		}
	}

	return this.button( this._groupSelector, idx );
} );

// Destroy the button sets selected
DataTable.Api.register( 'buttons().destroy()', function () {
	this.pluck( 'inst' ).unique().each( function ( inst ) {
		inst.destroy();
	} );

	return this;
} );

// Remove a button
DataTable.Api.registerPlural( 'buttons().remove()', 'buttons().remove()', function () {
	this.each( function ( set ) {
		set.inst.remove( set.node );
	} );

	return this;
} );

// Information box that can be used by buttons
var _infoTimer;
DataTable.Api.register( 'buttons.info()', function ( title, message, time ) {
	var that = this;

	if ( title === false ) {
		this.off('destroy.btn-info');
		_fadeOut(
			$('#datatables_buttons_info'),
			400,
			function () {
				$(this).remove();
			}
		);
		clearTimeout( _infoTimer );
		_infoTimer = null;

		return this;
	}

	if ( _infoTimer ) {
		clearTimeout( _infoTimer );
	}

	if ( $('#datatables_buttons_info').length ) {
		$('#datatables_buttons_info').remove();
	}

	title = title ? '<h2>'+title+'</h2>' : '';

	_fadeIn(
		$('<div id="datatables_buttons_info" class="dt-button-info"/>')
			.html( title )
			.append( $('<div/>')[ typeof message === 'string' ? 'html' : 'append' ]( message ) )
			.css( 'display', 'none' )
			.appendTo( 'body' )
	);

	if ( time !== undefined && time !== 0 ) {
		_infoTimer = setTimeout( function () {
			that.buttons.info( false );
		}, time );
	}

	this.on('destroy.btn-info', function () {
		that.buttons.info(false);
	});

	return this;
} );

// Get data from the table for export - this is common to a number of plug-in
// buttons so it is included in the Buttons core library
DataTable.Api.register( 'buttons.exportData()', function ( options ) {
	if ( this.context.length ) {
		return _exportData( new DataTable.Api( this.context[0] ), options );
	}
} );

// Get information about the export that is common to many of the export data
// types (DRY)
DataTable.Api.register( 'buttons.exportInfo()', function ( conf ) {
	if ( ! conf ) {
		conf = {};
	}

	return {
		filename: _filename( conf ),
		title: _title( conf ),
		messageTop: _message(this, conf.message || conf.messageTop, 'top'),
		messageBottom: _message(this, conf.messageBottom, 'bottom')
	};
} );



/**
 * Get the file name for an exported file.
 *
 * @param {object}	config Button configuration
 * @param {boolean} incExtension Include the file name extension
 */
var _filename = function ( config )
{
	// Backwards compatibility
	var filename = config.filename === '*' && config.title !== '*' && config.title !== undefined && config.title !== null && config.title !== '' ?
		config.title :
		config.filename;

	if ( typeof filename === 'function' ) {
		filename = filename();
	}

	if ( filename === undefined || filename === null ) {
		return null;
	}

	if ( filename.indexOf( '*' ) !== -1 ) {
		filename = filename.replace( '*', $('head > title').text() ).trim();
	}

	// Strip characters which the OS will object to
	filename = filename.replace(/[^a-zA-Z0-9_\u00A1-\uFFFF\.,\-_ !\(\)]/g, "");

	var extension = _stringOrFunction( config.extension );
	if ( ! extension ) {
		extension = '';
	}

	return filename + extension;
};

/**
 * Simply utility method to allow parameters to be given as a function
 *
 * @param {undefined|string|function} option Option
 * @return {null|string} Resolved value
 */
var _stringOrFunction = function ( option )
{
	if ( option === null || option === undefined ) {
		return null;
	}
	else if ( typeof option === 'function' ) {
		return option();
	}
	return option;
};

/**
 * Get the title for an exported file.
 *
 * @param {object} config	Button configuration
 */
var _title = function ( config )
{
	var title = _stringOrFunction( config.title );

	return title === null ?
		null : title.indexOf( '*' ) !== -1 ?
			title.replace( '*', $('head > title').text() || 'Exported data' ) :
			title;
};

var _message = function ( dt, option, position )
{
	var message = _stringOrFunction( option );
	if ( message === null ) {
		return null;
	}

	var caption = $('caption', dt.table().container()).eq(0);
	if ( message === '*' ) {
		var side = caption.css( 'caption-side' );
		if ( side !== position ) {
			return null;
		}

		return caption.length ?
			caption.text() :
			'';
	}

	return message;
};




var _exportTextarea = $('<textarea/>')[0];
var _exportData = function ( dt, inOpts )
{
	var config = $.extend( true, {}, {
		rows:           null,
		columns:        '',
		modifier:       {
			search: 'applied',
			order:  'applied'
		},
		orthogonal:     'display',
		stripHtml:      true,
		stripNewlines:  true,
		decodeEntities: true,
		trim:           true,
		format:         {
			header: function ( d ) {
				return Buttons.stripData( d, config );
			},
			footer: function ( d ) {
				return Buttons.stripData( d, config );
			},
			body: function ( d ) {
				return Buttons.stripData( d, config );
			}
		},
		customizeData: null
	}, inOpts );

	var header = dt.columns( config.columns ).indexes().map( function (idx) {
		var el = dt.column( idx ).header();
		return config.format.header( el.innerHTML, idx, el );
	} ).toArray();

	var footer = dt.table().footer() ?
		dt.columns( config.columns ).indexes().map( function (idx) {
			var el = dt.column( idx ).footer();
			return config.format.footer( el ? el.innerHTML : '', idx, el );
		} ).toArray() :
		null;
	
	// If Select is available on this table, and any rows are selected, limit the export
	// to the selected rows. If no rows are selected, all rows will be exported. Specify
	// a `selected` modifier to control directly.
	var modifier = $.extend( {}, config.modifier );
	if ( dt.select && typeof dt.select.info === 'function' && modifier.selected === undefined ) {
		if ( dt.rows( config.rows, $.extend( { selected: true }, modifier ) ).any() ) {
			$.extend( modifier, { selected: true } )
		}
	}

	var rowIndexes = dt.rows( config.rows, modifier ).indexes().toArray();
	var selectedCells = dt.cells( rowIndexes, config.columns );
	var cells = selectedCells
		.render( config.orthogonal )
		.toArray();
	var cellNodes = selectedCells
		.nodes()
		.toArray();

	var columns = header.length;
	var rows = columns > 0 ? cells.length / columns : 0;
	var body = [];
	var cellCounter = 0;

	for ( var i=0, ien=rows ; i<ien ; i++ ) {
		var row = [ columns ];

		for ( var j=0 ; j<columns ; j++ ) {
			row[j] = config.format.body( cells[ cellCounter ], i, j, cellNodes[ cellCounter ] );
			cellCounter++;
		}

		body[i] = row;
	}

	var data = {
		header: header,
		footer: footer,
		body:   body
	};

	if ( config.customizeData ) {
		config.customizeData( data );
	}

	return data;
};


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * DataTables interface
 */

// Attach to DataTables objects for global access
$.fn.dataTable.Buttons = Buttons;
$.fn.DataTable.Buttons = Buttons;



// DataTables creation - check if the buttons have been defined for this table,
// they will have been if the `B` option was used in `dom`, otherwise we should
// create the buttons instance here so they can be inserted into the document
// using the API. Listen for `init` for compatibility with pre 1.10.10, but to
// be removed in future.
$(document).on( 'init.dt plugin-init.dt', function (e, settings) {
	if ( e.namespace !== 'dt' ) {
		return;
	}

	var opts = settings.oInit.buttons || DataTable.defaults.buttons;

	if ( opts && ! settings._buttons ) {
		new Buttons( settings, opts ).container();
	}
} );

function _init ( settings, options ) {
	var api = new DataTable.Api( settings );
	var opts = options
		? options
		: api.init().buttons || DataTable.defaults.buttons;

	return new Buttons( api, opts ).container();
}

// DataTables `dom` feature option
DataTable.ext.feature.push( {
	fnInit: _init,
	cFeature: "B"
} );

// DataTables 2 layout feature
if ( DataTable.ext.features ) {
	DataTable.ext.features.register( 'buttons', _init );
}


return DataTable;
}));


/*! Bootstrap integration for DataTables' Buttons
 * ©2016 SpryMedia Ltd - datatables.net/license
 */

(function( factory ){
	if ( typeof define === 'function' && define.amd ) {
		// AMD
		define( ['jquery', 'datatables.net-bs5', 'datatables.net-buttons'], function ( $ ) {
			return factory( $, window, document );
		} );
	}
	else if ( typeof exports === 'object' ) {
		// CommonJS
		var jq = require('jquery');
		var cjsRequires = function (root, $) {
			if ( ! $.fn.dataTable ) {
				require('datatables.net-bs5')(root, $);
			}

			if ( ! $.fn.dataTable.Buttons ) {
				require('datatables.net-buttons')(root, $);
			}
		};

		if (typeof window !== 'undefined') {
			module.exports = function (root, $) {
				if ( ! root ) {
					// CommonJS environments without a window global must pass a
					// root. This will give an error otherwise
					root = window;
				}

				if ( ! $ ) {
					$ = jq( root );
				}

				cjsRequires( root, $ );
				return factory( $, root, root.document );
			};
		}
		else {
			cjsRequires( window, jq );
			module.exports = factory( jq, window, window.document );
		}
	}
	else {
		// Browser
		factory( jQuery, window, document );
	}
}(function( $, window, document, undefined ) {
'use strict';
var DataTable = $.fn.dataTable;



$.extend( true, DataTable.Buttons.defaults, {
	dom: {
		container: {
			className: 'dt-buttons btn-group flex-wrap'
		},
		button: {
			className: 'btn btn-secondary'
		},
		collection: {
			tag: 'div',
			className: 'dropdown-menu',
			closeButton: false,
			button: {
				tag: 'a',
				className: 'dt-button dropdown-item',
				active: 'active',
				disabled: 'disabled'
			}
		},
		splitWrapper: {
			tag: 'div',
			className: 'dt-btn-split-wrapper btn-group',
			closeButton: false,
		},
		splitDropdown: {
			tag: 'button',
			text: '',
			className: 'btn btn-secondary dt-btn-split-drop dropdown-toggle dropdown-toggle-split',
			closeButton: false,
			align: 'split-left',
			splitAlignClass: 'dt-button-split-left'
		},
		splitDropdownButton: {
			tag: 'button',
			className: 'dt-btn-split-drop-button btn btn-secondary',
			closeButton: false
		}
	},
	buttonCreated: function ( config, button ) {
		return config.buttons ?
			$('<div class="btn-group"/>').append(button) :
			button;
	}
} );

DataTable.ext.buttons.collection.className += ' dropdown-toggle';
DataTable.ext.buttons.collection.rightAlignClassName = 'dropdown-menu-right';


return DataTable;
}));


/*! SearchBuilder 1.4.2
 * ©SpryMedia Ltd - datatables.net/license/mit
 */

(function( factory ){
	if ( typeof define === 'function' && define.amd ) {
		// AMD
		define( ['jquery', 'datatables.net'], function ( $ ) {
			return factory( $, window, document );
		} );
	}
	else if ( typeof exports === 'object' ) {
		// CommonJS
		var jq = require('jquery');
		var cjsRequires = function (root, $) {
			if ( ! $.fn.dataTable ) {
				require('datatables.net')(root, $);
			}
		};

		if (typeof window !== 'undefined') {
			module.exports = function (root, $) {
				if ( ! root ) {
					// CommonJS environments without a window global must pass a
					// root. This will give an error otherwise
					root = window;
				}

				if ( ! $ ) {
					$ = jq( root );
				}

				cjsRequires( root, $ );
				return factory( $, root, root.document );
			};
		}
		else {
			cjsRequires( window, jq );
			module.exports = factory( jq, window, window.document );
		}
	}
	else {
		// Browser
		factory( jQuery, window, document );
	}
}(function( $, window, document, undefined ) {
'use strict';
var DataTable = $.fn.dataTable;


(function () {
    'use strict';

    var $$3;
    var dataTable$3;
    function moment() {
        return window.moment;
    }
    function luxon() {
        return window.luxon;
    }
    /**
     * Sets the value of jQuery for use in the file
     *
     * @param jq the instance of jQuery to be set
     */
    function setJQuery$2(jq) {
        $$3 = jq;
        dataTable$3 = jq.fn.dataTable;
    }
    /**
     * The Criteria class is used within SearchBuilder to represent a search criteria
     */
    var Criteria = /** @class */ (function () {
        function Criteria(table, opts, topGroup, index, depth, serverData) {
            if (index === void 0) { index = 0; }
            if (depth === void 0) { depth = 1; }
            if (serverData === void 0) { serverData = undefined; }
            var _this = this;
            // Check that the required version of DataTables is included
            if (!dataTable$3 || !dataTable$3.versionCheck || !dataTable$3.versionCheck('1.10.0')) {
                throw new Error('SearchPane requires DataTables 1.10 or newer');
            }
            this.classes = $$3.extend(true, {}, Criteria.classes);
            // Get options from user and any extra conditions/column types defined by plug-ins
            this.c = $$3.extend(true, {}, Criteria.defaults, $$3.fn.dataTable.ext.searchBuilder, opts);
            var i18n = this.c.i18n;
            this.s = {
                condition: undefined,
                conditions: {},
                data: undefined,
                dataIdx: -1,
                dataPoints: [],
                dateFormat: false,
                depth: depth,
                dt: table,
                filled: false,
                index: index,
                origData: undefined,
                preventRedraw: false,
                serverData: serverData,
                topGroup: topGroup,
                type: '',
                value: []
            };
            this.dom = {
                buttons: $$3('<div/>')
                    .addClass(this.classes.buttonContainer),
                condition: $$3('<select disabled/>')
                    .addClass(this.classes.condition)
                    .addClass(this.classes.dropDown)
                    .addClass(this.classes.italic)
                    .attr('autocomplete', 'hacking'),
                conditionTitle: $$3('<option value="" disabled selected hidden/>')
                    .html(this.s.dt.i18n('searchBuilder.condition', i18n.condition)),
                container: $$3('<div/>')
                    .addClass(this.classes.container),
                data: $$3('<select/>')
                    .addClass(this.classes.data)
                    .addClass(this.classes.dropDown)
                    .addClass(this.classes.italic),
                dataTitle: $$3('<option value="" disabled selected hidden/>')
                    .html(this.s.dt.i18n('searchBuilder.data', i18n.data)),
                defaultValue: $$3('<select disabled/>')
                    .addClass(this.classes.value)
                    .addClass(this.classes.dropDown)
                    .addClass(this.classes.select)
                    .addClass(this.classes.italic),
                "delete": $$3('<button/>')
                    .html(this.s.dt.i18n('searchBuilder.delete', i18n["delete"]))
                    .addClass(this.classes["delete"])
                    .addClass(this.classes.button)
                    .attr('title', this.s.dt.i18n('searchBuilder.deleteTitle', i18n.deleteTitle))
                    .attr('type', 'button'),
                inputCont: $$3('<div/>')
                    .addClass(this.classes.inputCont),
                // eslint-disable-next-line no-useless-escape
                left: $$3('<button/>')
                    .html(this.s.dt.i18n('searchBuilder.left', i18n.left))
                    .addClass(this.classes.left)
                    .addClass(this.classes.button)
                    .attr('title', this.s.dt.i18n('searchBuilder.leftTitle', i18n.leftTitle))
                    .attr('type', 'button'),
                // eslint-disable-next-line no-useless-escape
                right: $$3('<button/>')
                    .html(this.s.dt.i18n('searchBuilder.right', i18n.right))
                    .addClass(this.classes.right)
                    .addClass(this.classes.button)
                    .attr('title', this.s.dt.i18n('searchBuilder.rightTitle', i18n.rightTitle))
                    .attr('type', 'button'),
                value: [
                    $$3('<select disabled/>')
                        .addClass(this.classes.value)
                        .addClass(this.classes.dropDown)
                        .addClass(this.classes.italic)
                        .addClass(this.classes.select)
                ],
                valueTitle: $$3('<option value="--valueTitle--" disabled selected hidden/>')
                    .html(this.s.dt.i18n('searchBuilder.value', i18n.value))
            };
            // If the greyscale option is selected then add the class to add the grey colour to SearchBuilder
            if (this.c.greyscale) {
                this.dom.data.addClass(this.classes.greyscale);
                this.dom.condition.addClass(this.classes.greyscale);
                this.dom.defaultValue.addClass(this.classes.greyscale);
                for (var _i = 0, _a = this.dom.value; _i < _a.length; _i++) {
                    var val = _a[_i];
                    val.addClass(this.classes.greyscale);
                }
            }
            $$3(window).on('resize.dtsb', dataTable$3.util.throttle(function () {
                _this.s.topGroup.trigger('dtsb-redrawLogic');
            }));
            this._buildCriteria();
            return this;
        }
        /**
         * Escape html characters within a string
         *
         * @param txt the string to be escaped
         * @returns the escaped string
         */
        Criteria._escapeHTML = function (txt) {
            return txt
                .toString()
                .replace(/&amp;/g, '&')
                .replace(/&lt;/g, '<')
                .replace(/&gt;/g, '>')
                .replace(/&quot;/g, '"');
        };
        /**
         * Parses formatted numbers down to a form where they can be compared
         *
         * @param val the value to convert
         * @returns the converted value
         */
        Criteria.parseNumFmt = function (val) {
            return +val.replace(/(?!^-)[^0-9.]/g, '');
        };
        /**
         * Adds the left button to the criteria
         */
        Criteria.prototype.updateArrows = function (hasSiblings) {
            if (hasSiblings === void 0) { hasSiblings = false; }
            // Empty the container and append all of the elements in the correct order
            this.dom.container.children().detach();
            this.dom.container
                .append(this.dom.data)
                .append(this.dom.condition)
                .append(this.dom.inputCont);
            this.setListeners();
            // Trigger the inserted events for the value elements as they are inserted
            if (this.dom.value[0] !== undefined) {
                this.dom.value[0].trigger('dtsb-inserted');
            }
            for (var i = 1; i < this.dom.value.length; i++) {
                this.dom.inputCont.append(this.dom.value[i]);
                this.dom.value[i].trigger('dtsb-inserted');
            }
            // If this is a top level criteria then don't let it move left
            if (this.s.depth > 1) {
                this.dom.buttons.append(this.dom.left);
            }
            // If the depthLimit of the query has been hit then don't add the right button
            if ((this.c.depthLimit === false || this.s.depth < this.c.depthLimit) && hasSiblings) {
                this.dom.buttons.append(this.dom.right);
            }
            else {
                this.dom.right.remove();
            }
            this.dom.buttons.append(this.dom["delete"]);
            this.dom.container.append(this.dom.buttons);
        };
        /**
         * Destroys the criteria, removing listeners and container from the dom
         */
        Criteria.prototype.destroy = function () {
            // Turn off listeners
            this.dom.data.off('.dtsb');
            this.dom.condition.off('.dtsb');
            this.dom["delete"].off('.dtsb');
            for (var _i = 0, _a = this.dom.value; _i < _a.length; _i++) {
                var val = _a[_i];
                val.off('.dtsb');
            }
            // Remove container from the dom
            this.dom.container.remove();
        };
        /**
         * Passes in the data for the row and compares it against this single criteria
         *
         * @param rowData The data for the row to be compared
         * @returns boolean Whether the criteria has passed
         */
        Criteria.prototype.search = function (rowData, rowIdx) {
            var condition = this.s.conditions[this.s.condition];
            if (this.s.condition !== undefined && condition !== undefined) {
                var filter = rowData[this.s.dataIdx];
                // This check is in place for if a custom decimal character is in place
                if (this.s.type.includes('num') &&
                    (this.s.dt.settings()[0].oLanguage.sDecimal !== '' ||
                        this.s.dt.settings()[0].oLanguage.sThousands !== '')) {
                    var splitRD = [rowData[this.s.dataIdx]];
                    if (this.s.dt.settings()[0].oLanguage.sDecimal !== '') {
                        splitRD = rowData[this.s.dataIdx].split(this.s.dt.settings()[0].oLanguage.sDecimal);
                    }
                    if (this.s.dt.settings()[0].oLanguage.sThousands !== '') {
                        for (var i = 0; i < splitRD.length; i++) {
                            splitRD[i] = splitRD[i].replace(this.s.dt.settings()[0].oLanguage.sThousands, ',');
                        }
                    }
                    filter = splitRD.join('.');
                }
                // If orthogonal data is in place we need to get it's values for searching
                if (this.c.orthogonal.search !== 'filter') {
                    var settings = this.s.dt.settings()[0];
                    filter = settings.oApi._fnGetCellData(settings, rowIdx, this.s.dataIdx, typeof this.c.orthogonal === 'string' ?
                        this.c.orthogonal :
                        this.c.orthogonal.search);
                }
                if (this.s.type === 'array') {
                    // Make sure we are working with an array
                    if (!Array.isArray(filter)) {
                        filter = [filter];
                    }
                    filter.sort();
                    for (var _i = 0, filter_1 = filter; _i < filter_1.length; _i++) {
                        var filt = filter_1[_i];
                        if (filt && typeof filt === 'string') {
                            filt = filt.replace(/[\r\n\u2028]/g, ' ');
                        }
                    }
                }
                else if (filter !== null && typeof filter === 'string') {
                    filter = filter.replace(/[\r\n\u2028]/g, ' ');
                }
                if (this.s.type.includes('html') && typeof filter === 'string') {
                    filter = filter.replace(/(<([^>]+)>)/ig, '');
                }
                // Not ideal, but jqueries .val() returns an empty string even
                // when the value set is null, so we shall assume the two are equal
                if (filter === null) {
                    filter = '';
                }
                return condition.search(filter, this.s.value, this);
            }
        };
        /**
         * Gets the details required to rebuild the criteria
         */
        Criteria.prototype.getDetails = function (deFormatDates) {
            if (deFormatDates === void 0) { deFormatDates = false; }
            // This check is in place for if a custom decimal character is in place
            if (this.s.type !== null &&
                this.s.type.includes('num') &&
                (this.s.dt.settings()[0].oLanguage.sDecimal !== '' || this.s.dt.settings()[0].oLanguage.sThousands !== '')) {
                for (var i = 0; i < this.s.value.length; i++) {
                    var splitRD = [this.s.value[i].toString()];
                    if (this.s.dt.settings()[0].oLanguage.sDecimal !== '') {
                        splitRD = this.s.value[i].split(this.s.dt.settings()[0].oLanguage.sDecimal);
                    }
                    if (this.s.dt.settings()[0].oLanguage.sThousands !== '') {
                        for (var j = 0; j < splitRD.length; j++) {
                            splitRD[j] = splitRD[j].replace(this.s.dt.settings()[0].oLanguage.sThousands, ',');
                        }
                    }
                    this.s.value[i] = splitRD.join('.');
                }
            }
            else if (this.s.type !== null && deFormatDates) {
                if (this.s.type.includes('date') ||
                    this.s.type.includes('time')) {
                    for (var i = 0; i < this.s.value.length; i++) {
                        if (this.s.value[i].match(/^\d{4}-([0]\d|1[0-2])-([0-2]\d|3[01])$/g) === null) {
                            this.s.value[i] = '';
                        }
                    }
                }
                else if (this.s.type.includes('moment')) {
                    for (var i = 0; i < this.s.value.length; i++) {
                        if (this.s.value[i] &&
                            this.s.value[i].length > 0 &&
                            moment()(this.s.value[i], this.s.dateFormat, true).isValid()) {
                            this.s.value[i] = moment()(this.s.value[i], this.s.dateFormat).format('YYYY-MM-DD HH:mm:ss');
                        }
                    }
                }
                else if (this.s.type.includes('luxon')) {
                    for (var i = 0; i < this.s.value.length; i++) {
                        if (this.s.value[i] &&
                            this.s.value[i].length > 0 &&
                            luxon().DateTime.fromFormat(this.s.value[i], this.s.dateFormat).invalid === null) {
                            this.s.value[i] = luxon().DateTime.fromFormat(this.s.value[i], this.s.dateFormat).toFormat('yyyy-MM-dd HH:mm:ss');
                        }
                    }
                }
            }
            if (this.s.type.includes('num') && this.s.dt.page.info().serverSide) {
                for (var i = 0; i < this.s.value.length; i++) {
                    this.s.value[i] = this.s.value[i].replace(/[^0-9.\-]/g, '');
                }
            }
            return {
                condition: this.s.condition,
                data: this.s.data,
                origData: this.s.origData,
                type: this.s.type,
                value: this.s.value.map(function (a) { return a !== null && a !== undefined ? a.toString() : a; })
            };
        };
        /**
         * Getter for the node for the container of the criteria
         *
         * @returns JQuery<HTMLElement> the node for the container
         */
        Criteria.prototype.getNode = function () {
            return this.dom.container;
        };
        /**
         * Populates the criteria data, condition and value(s) as far as has been selected
         */
        Criteria.prototype.populate = function () {
            this._populateData();
            // If the column index has been found attempt to select a condition
            if (this.s.dataIdx !== -1) {
                this._populateCondition();
                // If the condittion has been found attempt to select the values
                if (this.s.condition !== undefined) {
                    this._populateValue();
                }
            }
        };
        /**
         * Rebuilds the criteria based upon the details passed in
         *
         * @param loadedCriteria the details required to rebuild the criteria
         */
        Criteria.prototype.rebuild = function (loadedCriteria) {
            // Check to see if the previously selected data exists, if so select it
            var foundData = false;
            var dataIdx;
            this._populateData();
            // If a data selection has previously been made attempt to find and select it
            if (loadedCriteria.data !== undefined) {
                var italic_1 = this.classes.italic;
                var data_1 = this.dom.data;
                this.dom.data.children('option').each(function () {
                    if (!foundData &&
                        ($$3(this).text() === loadedCriteria.data ||
                            loadedCriteria.origData && $$3(this).prop('origData') === loadedCriteria.origData)) {
                        $$3(this).prop('selected', true);
                        data_1.removeClass(italic_1);
                        foundData = true;
                        dataIdx = $$3(this).val();
                    }
                    else {
                        $$3(this).removeProp('selected');
                    }
                });
            }
            // If the data has been found and selected then the condition can be populated and searched
            if (foundData) {
                this.s.data = loadedCriteria.data;
                this.s.origData = loadedCriteria.origData;
                this.s.dataIdx = dataIdx;
                this.c.orthogonal = this._getOptions().orthogonal;
                this.dom.dataTitle.remove();
                this._populateCondition();
                this.dom.conditionTitle.remove();
                var condition = void 0;
                // Check to see if the previously selected condition exists, if so select it
                var options = this.dom.condition.children('option');
                // eslint-disable-next-line @typescript-eslint/prefer-for-of
                for (var i = 0; i < options.length; i++) {
                    var option = $$3(options[i]);
                    if (loadedCriteria.condition !== undefined &&
                        option.val() === loadedCriteria.condition &&
                        typeof loadedCriteria.condition === 'string') {
                        option.prop('selected', true);
                        condition = option.val();
                    }
                    else {
                        option.removeProp('selected');
                    }
                }
                this.s.condition = condition;
                // If the condition has been found and selected then the value can be populated and searched
                if (this.s.condition !== undefined) {
                    this.dom.conditionTitle.removeProp('selected');
                    this.dom.conditionTitle.remove();
                    this.dom.condition.removeClass(this.classes.italic);
                    // eslint-disable-next-line @typescript-eslint/prefer-for-of
                    for (var i = 0; i < options.length; i++) {
                        var option = $$3(options[i]);
                        if (option.val() !== this.s.condition) {
                            option.removeProp('selected');
                        }
                    }
                    this._populateValue(loadedCriteria);
                }
                else {
                    this.dom.conditionTitle.prependTo(this.dom.condition).prop('selected', true);
                }
            }
        };
        /**
         * Sets the listeners for the criteria
         */
        Criteria.prototype.setListeners = function () {
            var _this = this;
            this.dom.data
                .unbind('change')
                .on('change.dtsb', function () {
                _this.dom.dataTitle.removeProp('selected');
                // Need to go over every option to identify the correct selection
                var options = _this.dom.data.children('option.' + _this.classes.option);
                // eslint-disable-next-line @typescript-eslint/prefer-for-of
                for (var i = 0; i < options.length; i++) {
                    var option = $$3(options[i]);
                    if (option.val() === _this.dom.data.val()) {
                        _this.dom.data.removeClass(_this.classes.italic);
                        option.prop('selected', true);
                        _this.s.dataIdx = +option.val();
                        _this.s.data = option.text();
                        _this.s.origData = option.prop('origData');
                        _this.c.orthogonal = _this._getOptions().orthogonal;
                        // When the data is changed, the values in condition and
                        // value may also change so need to renew them
                        _this._clearCondition();
                        _this._clearValue();
                        _this._populateCondition();
                        // If this criteria was previously active in the search then
                        // remove it from the search and trigger a new search
                        if (_this.s.filled) {
                            _this.s.filled = false;
                            _this.s.dt.draw();
                            _this.setListeners();
                        }
                        _this.s.dt.state.save();
                    }
                    else {
                        option.removeProp('selected');
                    }
                }
            });
            this.dom.condition
                .unbind('change')
                .on('change.dtsb', function () {
                _this.dom.conditionTitle.removeProp('selected');
                // Need to go over every option to identify the correct selection
                var options = _this.dom.condition.children('option.' + _this.classes.option);
                // eslint-disable-next-line @typescript-eslint/prefer-for-of
                for (var i = 0; i < options.length; i++) {
                    var option = $$3(options[i]);
                    if (option.val() === _this.dom.condition.val()) {
                        _this.dom.condition.removeClass(_this.classes.italic);
                        option.prop('selected', true);
                        var condDisp = option.val();
                        // Find the condition that has been selected and store it internally
                        for (var _i = 0, _a = Object.keys(_this.s.conditions); _i < _a.length; _i++) {
                            var cond = _a[_i];
                            if (cond === condDisp) {
                                _this.s.condition = condDisp;
                                break;
                            }
                        }
                        // When the condition is changed, the value selector may switch between
                        // a select element and an input element
                        _this._clearValue();
                        _this._populateValue();
                        for (var _b = 0, _c = _this.dom.value; _b < _c.length; _b++) {
                            var val = _c[_b];
                            // If this criteria was previously active in the search then remove
                            // it from the search and trigger a new search
                            if (_this.s.filled && val !== undefined && _this.dom.inputCont.has(val[0]).length !== 0) {
                                _this.s.filled = false;
                                _this.s.dt.draw();
                                _this.setListeners();
                            }
                        }
                        if (_this.dom.value.length === 0 ||
                            _this.dom.value.length === 1 && _this.dom.value[0] === undefined) {
                            _this.s.dt.draw();
                        }
                    }
                    else {
                        option.removeProp('selected');
                    }
                }
            });
        };
        Criteria.prototype.setupButtons = function () {
            if (window.innerWidth > 550) {
                this.dom.container.removeClass(this.classes.vertical);
                this.dom.buttons.css('left', null);
                this.dom.buttons.css('top', null);
                return;
            }
            this.dom.container.addClass(this.classes.vertical);
            this.dom.buttons.css('left', this.dom.data.innerWidth());
            this.dom.buttons.css('top', this.dom.data.position().top);
        };
        /**
         * Builds the elements of the dom together
         */
        Criteria.prototype._buildCriteria = function () {
            // Append Titles for select elements
            this.dom.data.append(this.dom.dataTitle);
            this.dom.condition.append(this.dom.conditionTitle);
            // Add elements to container
            this.dom.container
                .append(this.dom.data)
                .append(this.dom.condition);
            this.dom.inputCont.empty();
            for (var _i = 0, _a = this.dom.value; _i < _a.length; _i++) {
                var val = _a[_i];
                val.append(this.dom.valueTitle);
                this.dom.inputCont.append(val);
            }
            // Add buttons to container
            this.dom.buttons
                .append(this.dom["delete"])
                .append(this.dom.right);
            this.dom.container.append(this.dom.inputCont).append(this.dom.buttons);
            this.setListeners();
        };
        /**
         * Clears the condition select element
         */
        Criteria.prototype._clearCondition = function () {
            this.dom.condition.empty();
            this.dom.conditionTitle.prop('selected', true).attr('disabled', 'true');
            this.dom.condition.prepend(this.dom.conditionTitle).prop('selectedIndex', 0);
            this.s.conditions = {};
            this.s.condition = undefined;
        };
        /**
         * Clears the value elements
         */
        Criteria.prototype._clearValue = function () {
            if (this.s.condition !== undefined) {
                if (this.dom.value.length > 0 && this.dom.value[0] !== undefined) {
                    var _loop_1 = function (val) {
                        if (val !== undefined) {
                            // Timeout is annoying but because of IOS
                            setTimeout(function () {
                                val.remove();
                            }, 50);
                        }
                    };
                    // Remove all of the value elements
                    for (var _i = 0, _a = this.dom.value; _i < _a.length; _i++) {
                        var val = _a[_i];
                        _loop_1(val);
                    }
                }
                // Call the init function to get the value elements for this condition
                this.dom.value = [].concat(this.s.conditions[this.s.condition].init(this, Criteria.updateListener));
                if (this.dom.value.length > 0 && this.dom.value[0] !== undefined) {
                    this.dom.inputCont
                        .empty()
                        .append(this.dom.value[0])
                        .insertAfter(this.dom.condition);
                    this.dom.value[0].trigger('dtsb-inserted');
                    // Insert all of the value elements
                    for (var i = 1; i < this.dom.value.length; i++) {
                        this.dom.inputCont.append(this.dom.value[i]);
                        this.dom.value[i].trigger('dtsb-inserted');
                    }
                }
            }
            else {
                var _loop_2 = function (val) {
                    if (val !== undefined) {
                        // Timeout is annoying but because of IOS
                        setTimeout(function () {
                            val.remove();
                        }, 50);
                    }
                };
                // Remove all of the value elements
                for (var _b = 0, _c = this.dom.value; _b < _c.length; _b++) {
                    var val = _c[_b];
                    _loop_2(val);
                }
                // Append the default valueTitle to the default select element
                this.dom.valueTitle
                    .prop('selected', true);
                this.dom.defaultValue
                    .append(this.dom.valueTitle)
                    .insertAfter(this.dom.condition);
            }
            this.s.value = [];
            this.dom.value = [
                $$3('<select disabled/>')
                    .addClass(this.classes.value)
                    .addClass(this.classes.dropDown)
                    .addClass(this.classes.italic)
                    .addClass(this.classes.select)
                    .append(this.dom.valueTitle.clone())
            ];
        };
        /**
         * Gets the options for the column
         *
         * @returns {object} The options for the column
         */
        Criteria.prototype._getOptions = function () {
            var table = this.s.dt;
            return $$3.extend(true, {}, Criteria.defaults, table.settings()[0].aoColumns[this.s.dataIdx].searchBuilder);
        };
        /**
         * Populates the condition dropdown
         */
        Criteria.prototype._populateCondition = function () {
            var conditionOpts = [];
            var conditionsLength = Object.keys(this.s.conditions).length;
            var colInits = this.s.dt.settings()[0].aoColumns;
            var column = +this.dom.data.children('option:selected').val();
            // If there are no conditions stored then we need to get them from the appropriate type
            if (conditionsLength === 0) {
                this.s.type = this.s.dt.columns().type().toArray()[column];
                if (colInits !== undefined) {
                    var colInit = colInits[column];
                    if (colInit.searchBuilderType !== undefined && colInit.searchBuilderType !== null) {
                        this.s.type = colInit.searchBuilderType;
                    }
                    else if (this.s.type === undefined || this.s.type === null) {
                        this.s.type = colInit.sType;
                    }
                }
                // If the column type is still unknown, call a draw to try reading it again
                if (this.s.type === null || this.s.type === undefined) {
                    $$3.fn.dataTable.ext.oApi._fnColumnTypes(this.s.dt.settings()[0]);
                    this.s.type = this.s.dt.columns().type().toArray()[column];
                }
                // Enable the condition element
                this.dom.condition
                    .removeAttr('disabled')
                    .empty()
                    .append(this.dom.conditionTitle)
                    .addClass(this.classes.italic);
                this.dom.conditionTitle
                    .prop('selected', true);
                var decimal = this.s.dt.settings()[0].oLanguage.sDecimal;
                // This check is in place for if a custom decimal character is in place
                if (decimal !== '' && this.s.type.indexOf(decimal) === this.s.type.length - decimal.length) {
                    if (this.s.type.includes('num-fmt')) {
                        this.s.type = this.s.type.replace(decimal, '');
                    }
                    else if (this.s.type.includes('num')) {
                        this.s.type = this.s.type.replace(decimal, '');
                    }
                }
                // Select which conditions are going to be used based on the column type
                var conditionObj = this.c.conditions[this.s.type] !== undefined ?
                    this.c.conditions[this.s.type] :
                    this.s.type.includes('moment') ?
                        this.c.conditions.moment :
                        this.s.type.includes('luxon') ?
                            this.c.conditions.luxon :
                            this.c.conditions.string;
                // If it is a moment format then extract the date format
                if (this.s.type.includes('moment')) {
                    this.s.dateFormat = this.s.type.replace(/moment-/g, '');
                }
                else if (this.s.type.includes('luxon')) {
                    this.s.dateFormat = this.s.type.replace(/luxon-/g, '');
                }
                // Add all of the conditions to the select element
                for (var _i = 0, _a = Object.keys(conditionObj); _i < _a.length; _i++) {
                    var condition = _a[_i];
                    if (conditionObj[condition] !== null) {
                        // Serverside processing does not supply the options for the select elements
                        // Instead input elements need to be used for these instead
                        if (this.s.dt.page.info().serverSide && conditionObj[condition].init === Criteria.initSelect) {
                            var col = colInits[column];
                            if (this.s.serverData && this.s.serverData[col.data]) {
                                conditionObj[condition].init = Criteria.initSelectSSP;
                                conditionObj[condition].inputValue = Criteria.inputValueSelect;
                                conditionObj[condition].isInputValid = Criteria.isInputValidSelect;
                            }
                            else {
                                conditionObj[condition].init = Criteria.initInput;
                                conditionObj[condition].inputValue = Criteria.inputValueInput;
                                conditionObj[condition].isInputValid = Criteria.isInputValidInput;
                            }
                        }
                        this.s.conditions[condition] = conditionObj[condition];
                        var condName = conditionObj[condition].conditionName;
                        if (typeof condName === 'function') {
                            condName = condName(this.s.dt, this.c.i18n);
                        }
                        conditionOpts.push($$3('<option>', {
                            text: condName,
                            value: condition
                        })
                            .addClass(this.classes.option)
                            .addClass(this.classes.notItalic));
                    }
                }
            }
            // Otherwise we can just load them in
            else if (conditionsLength > 0) {
                this.dom.condition.empty().removeAttr('disabled').addClass(this.classes.italic);
                for (var _b = 0, _c = Object.keys(this.s.conditions); _b < _c.length; _b++) {
                    var condition = _c[_b];
                    var condName = this.s.conditions[condition].conditionName;
                    if (typeof condName === 'function') {
                        condName = condName(this.s.dt, this.c.i18n);
                    }
                    var newOpt = $$3('<option>', {
                        text: condName,
                        value: condition
                    })
                        .addClass(this.classes.option)
                        .addClass(this.classes.notItalic);
                    if (this.s.condition !== undefined && this.s.condition === condName) {
                        newOpt.prop('selected', true);
                        this.dom.condition.removeClass(this.classes.italic);
                    }
                    conditionOpts.push(newOpt);
                }
            }
            else {
                this.dom.condition
                    .attr('disabled', 'true')
                    .addClass(this.classes.italic);
                return;
            }
            for (var _d = 0, conditionOpts_1 = conditionOpts; _d < conditionOpts_1.length; _d++) {
                var opt = conditionOpts_1[_d];
                this.dom.condition.append(opt);
            }
            // Selecting a default condition if one is set
            if (colInits[column].searchBuilder && colInits[column].searchBuilder.defaultCondition) {
                var defaultCondition = colInits[column].searchBuilder.defaultCondition;
                // If it is a number just use it as an index
                if (typeof defaultCondition === 'number') {
                    this.dom.condition.prop('selectedIndex', defaultCondition);
                    this.dom.condition.trigger('change');
                }
                // If it is a string then things get slightly more tricly
                else if (typeof defaultCondition === 'string') {
                    // We need to check each condition option to see if any will match
                    for (var i = 0; i < conditionOpts.length; i++) {
                        // Need to check against the stored conditions so we can match the token "cond" to the option
                        for (var _e = 0, _f = Object.keys(this.s.conditions); _e < _f.length; _e++) {
                            var cond = _f[_e];
                            var condName = this.s.conditions[cond].conditionName;
                            if (
                            // If the conditionName matches the text of the option
                            (typeof condName === 'string' ? condName : condName(this.s.dt, this.c.i18n)) ===
                                conditionOpts[i].text() &&
                                // and the tokens match
                                cond === defaultCondition) {
                                // Select that option
                                this.dom.condition
                                    .prop('selectedIndex', this.dom.condition.children().toArray().indexOf(conditionOpts[i][0]))
                                    .removeClass(this.classes.italic);
                                this.dom.condition.trigger('change');
                                i = conditionOpts.length;
                                break;
                            }
                        }
                    }
                }
            }
            // If not default set then default to 0, the title
            else {
                this.dom.condition.prop('selectedIndex', 0);
            }
        };
        /**
         * Populates the data / column select element
         */
        Criteria.prototype._populateData = function () {
            var columns = this.s.dt.settings()[0].aoColumns;
            var includeColumns = this.s.dt.columns(this.c.columns).indexes().toArray();
            this.dom.data.empty().append(this.dom.dataTitle);
            for (var index = 0; index < columns.length; index++) {
                // Need to check that the column can be filtered on before adding it
                if (this.c.columns === true || includeColumns.includes(index)) {
                    var col = columns[index];
                    var opt = {
                        index: index,
                        origData: col.data,
                        text: (col.searchBuilderTitle || col.sTitle)
                            .replace(/(<([^>]+)>)/ig, '')
                    };
                    this.dom.data.append($$3('<option>', {
                        text: opt.text,
                        value: opt.index
                    })
                        .addClass(this.classes.option)
                        .addClass(this.classes.notItalic)
                        .prop('origData', col.data)
                        .prop('selected', this.s.dataIdx === opt.index ? true : false));
                    if (this.s.dataIdx === opt.index) {
                        this.dom.dataTitle.removeProp('selected');
                    }
                }
            }
        };
        /**
         * Populates the Value select element
         *
         * @param loadedCriteria optional, used to reload criteria from predefined filters
         */
        Criteria.prototype._populateValue = function (loadedCriteria) {
            var _this = this;
            var prevFilled = this.s.filled;
            this.s.filled = false;
            // Remove any previous value elements
            // Timeout is annoying but because of IOS
            setTimeout(function () {
                _this.dom.defaultValue.remove();
            }, 50);
            var _loop_3 = function (val) {
                // Timeout is annoying but because of IOS
                setTimeout(function () {
                    if (val !== undefined) {
                        val.remove();
                    }
                }, 50);
            };
            for (var _i = 0, _a = this.dom.value; _i < _a.length; _i++) {
                var val = _a[_i];
                _loop_3(val);
            }
            var children = this.dom.inputCont.children();
            if (children.length > 1) {
                // eslint-disable-next-line @typescript-eslint/prefer-for-of
                for (var i = 0; i < children.length; i++) {
                    $$3(children[i]).remove();
                }
            }
            // Find the column with the title matching the data for the criteria and take note of the index
            if (loadedCriteria !== undefined) {
                this.s.dt.columns().every(function (index) {
                    if (_this.s.dt.settings()[0].aoColumns[index].sTitle === loadedCriteria.data) {
                        _this.s.dataIdx = index;
                    }
                });
            }
            // Initialise the value elements based on the condition
            this.dom.value = [].concat(this.s.conditions[this.s.condition].init(this, Criteria.updateListener, loadedCriteria !== undefined ? loadedCriteria.value : undefined));
            if (loadedCriteria !== undefined && loadedCriteria.value !== undefined) {
                this.s.value = loadedCriteria.value;
            }
            this.dom.inputCont.empty();
            // Insert value elements and trigger the inserted event
            if (this.dom.value[0] !== undefined) {
                this.dom.value[0]
                    .appendTo(this.dom.inputCont)
                    .trigger('dtsb-inserted');
            }
            for (var i = 1; i < this.dom.value.length; i++) {
                this.dom.value[i]
                    .insertAfter(this.dom.value[i - 1])
                    .trigger('dtsb-inserted');
            }
            // Check if the criteria can be used in a search
            this.s.filled = this.s.conditions[this.s.condition].isInputValid(this.dom.value, this);
            this.setListeners();
            // If it can and this is different to before then trigger a draw
            if (!this.s.preventRedraw && prevFilled !== this.s.filled) {
                // If using SSP we want to restrict the amount of server calls that take place
                //  and this will already have taken place
                if (!this.s.dt.page.info().serverSide) {
                    this.s.dt.draw();
                }
                this.setListeners();
            }
        };
        /**
         * Provides throttling capabilities to SearchBuilder without having to use dt's _fnThrottle function
         * This is because that function is not quite suitable for our needs as it runs initially rather than waiting
         *
         * @param args arguments supplied to the throttle function
         * @returns Function that is to be run that implements the throttling
         */
        Criteria.prototype._throttle = function (fn, frequency) {
            if (frequency === void 0) { frequency = 200; }
            var last = null;
            var timer = null;
            var that = this;
            if (frequency === null) {
                frequency = 200;
            }
            return function () {
                var args = [];
                for (var _i = 0; _i < arguments.length; _i++) {
                    args[_i] = arguments[_i];
                }
                var now = +new Date();
                if (last !== null && now < last + frequency) {
                    clearTimeout(timer);
                }
                else {
                    last = now;
                }
                timer = setTimeout(function () {
                    last = null;
                    fn.apply(that, args);
                }, frequency);
            };
        };
        Criteria.version = '1.1.0';
        Criteria.classes = {
            button: 'dtsb-button',
            buttonContainer: 'dtsb-buttonContainer',
            condition: 'dtsb-condition',
            container: 'dtsb-criteria',
            data: 'dtsb-data',
            "delete": 'dtsb-delete',
            dropDown: 'dtsb-dropDown',
            greyscale: 'dtsb-greyscale',
            input: 'dtsb-input',
            inputCont: 'dtsb-inputCont',
            italic: 'dtsb-italic',
            joiner: 'dtsp-joiner',
            left: 'dtsb-left',
            notItalic: 'dtsb-notItalic',
            option: 'dtsb-option',
            right: 'dtsb-right',
            select: 'dtsb-select',
            value: 'dtsb-value',
            vertical: 'dtsb-vertical'
        };
        /**
         * Default initialisation function for select conditions
         */
        Criteria.initSelect = function (that, fn, preDefined, array) {
            if (preDefined === void 0) { preDefined = null; }
            if (array === void 0) { array = false; }
            var column = that.dom.data.children('option:selected').val();
            var indexArray = that.s.dt.rows().indexes().toArray();
            var settings = that.s.dt.settings()[0];
            that.dom.valueTitle.prop('selected', true);
            // Declare select element to be used with all of the default classes and listeners.
            var el = $$3('<select/>')
                .addClass(Criteria.classes.value)
                .addClass(Criteria.classes.dropDown)
                .addClass(Criteria.classes.italic)
                .addClass(Criteria.classes.select)
                .append(that.dom.valueTitle)
                .on('change.dtsb', function () {
                $$3(this).removeClass(Criteria.classes.italic);
                fn(that, this);
            });
            if (that.c.greyscale) {
                el.addClass(Criteria.classes.greyscale);
            }
            var added = [];
            var options = [];
            // Add all of the options from the table to the select element.
            // Only add one option for each possible value
            for (var _i = 0, indexArray_1 = indexArray; _i < indexArray_1.length; _i++) {
                var index = indexArray_1[_i];
                var filter = settings.oApi._fnGetCellData(settings, index, column, typeof that.c.orthogonal === 'string' ?
                    that.c.orthogonal :
                    that.c.orthogonal.search);
                var value = {
                    filter: typeof filter === 'string' ?
                        filter.replace(/[\r\n\u2028]/g, ' ') : // Need to replace certain characters to match search values
                        filter,
                    index: index,
                    text: settings.oApi._fnGetCellData(settings, index, column, typeof that.c.orthogonal === 'string' ?
                        that.c.orthogonal :
                        that.c.orthogonal.display)
                };
                // If we are dealing with an array type, either make sure we are working with arrays, or sort them
                if (that.s.type === 'array') {
                    value.filter = !Array.isArray(value.filter) ? [value.filter] : value.filter;
                    value.text = !Array.isArray(value.text) ? [value.text] : value.text;
                }
                // Function to add an option to the select element
                var addOption = function (filt, text) {
                    if (that.s.type.includes('html') && filt !== null && typeof filt === 'string') {
                        filt.replace(/(<([^>]+)>)/ig, '');
                    }
                    // Add text and value, stripping out any html if that is the column type
                    var opt = $$3('<option>', {
                        type: Array.isArray(filt) ? 'Array' : 'String',
                        value: filt
                    })
                        .data('sbv', filt)
                        .addClass(that.classes.option)
                        .addClass(that.classes.notItalic)
                        // Have to add the text this way so that special html characters are not escaped - &amp; etc.
                        .html(typeof text === 'string' ?
                        text.replace(/(<([^>]+)>)/ig, '') :
                        text);
                    var val = opt.val();
                    // Check that this value has not already been added
                    if (added.indexOf(val) === -1) {
                        added.push(val);
                        options.push(opt);
                        if (preDefined !== null && Array.isArray(preDefined[0])) {
                            preDefined[0] = preDefined[0].sort().join(',');
                        }
                        // If this value was previously selected as indicated by preDefined, then select it again
                        if (preDefined !== null && opt.val() === preDefined[0]) {
                            opt.prop('selected', true);
                            el.removeClass(Criteria.classes.italic);
                            that.dom.valueTitle.removeProp('selected');
                        }
                    }
                };
                // If this is to add the individual values within the array we need to loop over the array
                if (array) {
                    for (var i = 0; i < value.filter.length; i++) {
                        addOption(value.filter[i], value.text[i]);
                    }
                }
                // Otherwise the value that is in the cell is to be added
                else {
                    addOption(value.filter, Array.isArray(value.text) ? value.text.join(', ') : value.text);
                }
            }
            options.sort(function (a, b) {
                if (that.s.type === 'array' ||
                    that.s.type === 'string' ||
                    that.s.type === 'html') {
                    if (a.val() < b.val()) {
                        return -1;
                    }
                    else if (a.val() > b.val()) {
                        return 1;
                    }
                    else {
                        return 0;
                    }
                }
                else if (that.s.type === 'num' ||
                    that.s.type === 'html-num') {
                    if (+a.val().replace(/(<([^>]+)>)/ig, '') < +b.val().replace(/(<([^>]+)>)/ig, '')) {
                        return -1;
                    }
                    else if (+a.val().replace(/(<([^>]+)>)/ig, '') > +b.val().replace(/(<([^>]+)>)/ig, '')) {
                        return 1;
                    }
                    else {
                        return 0;
                    }
                }
                else if (that.s.type === 'num-fmt' || that.s.type === 'html-num-fmt') {
                    if (+a.val().replace(/[^0-9.]/g, '') < +b.val().replace(/[^0-9.]/g, '')) {
                        return -1;
                    }
                    else if (+a.val().replace(/[^0-9.]/g, '') > +b.val().replace(/[^0-9.]/g, '')) {
                        return 1;
                    }
                    else {
                        return 0;
                    }
                }
            });
            for (var _a = 0, options_1 = options; _a < options_1.length; _a++) {
                var opt = options_1[_a];
                el.append(opt);
            }
            return el;
        };
        /**
         * Default initialisation function for select conditions
         */
        Criteria.initSelectSSP = function (that, fn, preDefined) {
            if (preDefined === void 0) { preDefined = null; }
            that.dom.valueTitle.prop('selected', true);
            // Declare select element to be used with all of the default classes and listeners.
            var el = $$3('<select/>')
                .addClass(Criteria.classes.value)
                .addClass(Criteria.classes.dropDown)
                .addClass(Criteria.classes.italic)
                .addClass(Criteria.classes.select)
                .append(that.dom.valueTitle)
                .on('change.dtsb', function () {
                $$3(this).removeClass(Criteria.classes.italic);
                fn(that, this);
            });
            if (that.c.greyscale) {
                el.addClass(Criteria.classes.greyscale);
            }
            var options = [];
            for (var _i = 0, _a = that.s.serverData[that.s.origData]; _i < _a.length; _i++) {
                var option = _a[_i];
                var value = option.value;
                var label = option.label;
                // Function to add an option to the select element
                var addOption = function (filt, text) {
                    if (that.s.type.includes('html') && filt !== null && typeof filt === 'string') {
                        filt.replace(/(<([^>]+)>)/ig, '');
                    }
                    // Add text and value, stripping out any html if that is the column type
                    var opt = $$3('<option>', {
                        type: Array.isArray(filt) ? 'Array' : 'String',
                        value: filt
                    })
                        .data('sbv', filt)
                        .addClass(that.classes.option)
                        .addClass(that.classes.notItalic)
                        // Have to add the text this way so that special html characters are not escaped - &amp; etc.
                        .html(typeof text === 'string' ?
                        text.replace(/(<([^>]+)>)/ig, '') :
                        text);
                    options.push(opt);
                    // If this value was previously selected as indicated by preDefined, then select it again
                    if (preDefined !== null && opt.val() === preDefined[0]) {
                        opt.prop('selected', true);
                        el.removeClass(Criteria.classes.italic);
                        that.dom.valueTitle.removeProp('selected');
                    }
                };
                addOption(value, label);
            }
            for (var _b = 0, options_2 = options; _b < options_2.length; _b++) {
                var opt = options_2[_b];
                el.append(opt);
            }
            return el;
        };
        /**
         * Default initialisation function for select array conditions
         *
         * This exists because there needs to be different select functionality for contains/without and equals/not
         */
        Criteria.initSelectArray = function (that, fn, preDefined) {
            if (preDefined === void 0) { preDefined = null; }
            return Criteria.initSelect(that, fn, preDefined, true);
        };
        /**
         * Default initialisation function for input conditions
         */
        Criteria.initInput = function (that, fn, preDefined) {
            if (preDefined === void 0) { preDefined = null; }
            // Declare the input element
            var searchDelay = that.s.dt.settings()[0].searchDelay;
            var el = $$3('<input/>')
                .addClass(Criteria.classes.value)
                .addClass(Criteria.classes.input)
                .on('input.dtsb keypress.dtsb', that._throttle(function (e) {
                var code = e.keyCode || e.which;
                return fn(that, this, code);
            }, searchDelay === null ? 100 : searchDelay));
            if (that.c.greyscale) {
                el.addClass(Criteria.classes.greyscale);
            }
            // If there is a preDefined value then add it
            if (preDefined !== null) {
                el.val(preDefined[0]);
            }
            // This is add responsive functionality to the logic button without redrawing everything else
            that.s.dt.one('draw.dtsb', function () {
                that.s.topGroup.trigger('dtsb-redrawLogic');
            });
            return el;
        };
        /**
         * Default initialisation function for conditions requiring 2 inputs
         */
        Criteria.init2Input = function (that, fn, preDefined) {
            if (preDefined === void 0) { preDefined = null; }
            // Declare all of the necessary jQuery elements
            var searchDelay = that.s.dt.settings()[0].searchDelay;
            var els = [
                $$3('<input/>')
                    .addClass(Criteria.classes.value)
                    .addClass(Criteria.classes.input)
                    .on('input.dtsb keypress.dtsb', that._throttle(function (e) {
                    var code = e.keyCode || e.which;
                    return fn(that, this, code);
                }, searchDelay === null ? 100 : searchDelay)),
                $$3('<span>')
                    .addClass(that.classes.joiner)
                    .html(that.s.dt.i18n('searchBuilder.valueJoiner', that.c.i18n.valueJoiner)),
                $$3('<input/>')
                    .addClass(Criteria.classes.value)
                    .addClass(Criteria.classes.input)
                    .on('input.dtsb keypress.dtsb', that._throttle(function (e) {
                    var code = e.keyCode || e.which;
                    return fn(that, this, code);
                }, searchDelay === null ? 100 : searchDelay))
            ];
            if (that.c.greyscale) {
                els[0].addClass(Criteria.classes.greyscale);
                els[2].addClass(Criteria.classes.greyscale);
            }
            // If there is a preDefined value then add it
            if (preDefined !== null) {
                els[0].val(preDefined[0]);
                els[2].val(preDefined[1]);
            }
            // This is add responsive functionality to the logic button without redrawing everything else
            that.s.dt.one('draw.dtsb', function () {
                that.s.topGroup.trigger('dtsb-redrawLogic');
            });
            return els;
        };
        /**
         * Default initialisation function for date conditions
         */
        Criteria.initDate = function (that, fn, preDefined) {
            if (preDefined === void 0) { preDefined = null; }
            var searchDelay = that.s.dt.settings()[0].searchDelay;
            // Declare date element using DataTables dateTime plugin
            var el = $$3('<input/>')
                .addClass(Criteria.classes.value)
                .addClass(Criteria.classes.input)
                .dtDateTime({
                attachTo: 'input',
                format: that.s.dateFormat ? that.s.dateFormat : undefined
            })
                .on('change.dtsb', that._throttle(function () {
                return fn(that, this);
            }, searchDelay === null ? 100 : searchDelay))
                .on('input.dtsb keypress.dtsb', function (e) {
                that._throttle(function () {
                    var code = e.keyCode || e.which;
                    return fn(that, this, code);
                }, searchDelay === null ? 100 : searchDelay);
            });
            if (that.c.greyscale) {
                el.addClass(Criteria.classes.greyscale);
            }
            // If there is a preDefined value then add it
            if (preDefined !== null) {
                el.val(preDefined[0]);
            }
            // This is add responsive functionality to the logic button without redrawing everything else
            that.s.dt.one('draw.dtsb', function () {
                that.s.topGroup.trigger('dtsb-redrawLogic');
            });
            return el;
        };
        Criteria.initNoValue = function (that) {
            // This is add responsive functionality to the logic button without redrawing everything else
            that.s.dt.one('draw.dtsb', function () {
                that.s.topGroup.trigger('dtsb-redrawLogic');
            });
            return [];
        };
        Criteria.init2Date = function (that, fn, preDefined) {
            var _this = this;
            if (preDefined === void 0) { preDefined = null; }
            var searchDelay = that.s.dt.settings()[0].searchDelay;
            // Declare all of the date elements that are required using DataTables dateTime plugin
            var els = [
                $$3('<input/>')
                    .addClass(Criteria.classes.value)
                    .addClass(Criteria.classes.input)
                    .dtDateTime({
                    attachTo: 'input',
                    format: that.s.dateFormat ? that.s.dateFormat : undefined
                })
                    .on('change.dtsb', searchDelay !== null ?
                    that.s.dt.settings()[0].oApi._fnThrottle(function () {
                        return fn(that, this);
                    }, searchDelay) :
                    function () {
                        fn(that, _this);
                    })
                    .on('input.dtsb keypress.dtsb', function (e) {
                    that.s.dt.settings()[0].oApi._fnThrottle(function () {
                        var code = e.keyCode || e.which;
                        return fn(that, this, code);
                    }, searchDelay === null ? 0 : searchDelay);
                }),
                $$3('<span>')
                    .addClass(that.classes.joiner)
                    .html(that.s.dt.i18n('searchBuilder.valueJoiner', that.c.i18n.valueJoiner)),
                $$3('<input/>')
                    .addClass(Criteria.classes.value)
                    .addClass(Criteria.classes.input)
                    .dtDateTime({
                    attachTo: 'input',
                    format: that.s.dateFormat ? that.s.dateFormat : undefined
                })
                    .on('change.dtsb', searchDelay !== null ?
                    that.s.dt.settings()[0].oApi._fnThrottle(function () {
                        return fn(that, this);
                    }, searchDelay) :
                    function () {
                        fn(that, _this);
                    })
                    .on('input.dtsb keypress.dtsb', !that.c.enterSearch &&
                    !(that.s.dt.settings()[0].oInit.search !== undefined &&
                        that.s.dt.settings()[0].oInit.search["return"]) &&
                    searchDelay !== null ?
                    that.s.dt.settings()[0].oApi._fnThrottle(function () {
                        return fn(that, this);
                    }, searchDelay) :
                    function (e) {
                        var code = e.keyCode || e.which;
                        fn(that, _this, code);
                    })
            ];
            if (that.c.greyscale) {
                els[0].addClass(Criteria.classes.greyscale);
                els[2].addClass(Criteria.classes.greyscale);
            }
            // If there are and preDefined values then add them
            if (preDefined !== null && preDefined.length > 0) {
                els[0].val(preDefined[0]);
                els[2].val(preDefined[1]);
            }
            // This is add responsive functionality to the logic button without redrawing everything else
            that.s.dt.one('draw.dtsb', function () {
                that.s.topGroup.trigger('dtsb-redrawLogic');
            });
            return els;
        };
        /**
         * Default function for select elements to validate condition
         */
        Criteria.isInputValidSelect = function (el) {
            var allFilled = true;
            // Check each element to make sure that the selections are valid
            for (var _i = 0, el_1 = el; _i < el_1.length; _i++) {
                var element = el_1[_i];
                if (element.children('option:selected').length ===
                    element.children('option').length -
                        element.children('option.' + Criteria.classes.notItalic).length &&
                    element.children('option:selected').length === 1 &&
                    element.children('option:selected')[0] === element.children('option')[0]) {
                    allFilled = false;
                }
            }
            return allFilled;
        };
        /**
         * Default function for input and date elements to validate condition
         */
        Criteria.isInputValidInput = function (el) {
            var allFilled = true;
            // Check each element to make sure that the inputs are valid
            for (var _i = 0, el_2 = el; _i < el_2.length; _i++) {
                var element = el_2[_i];
                if (element.is('input') && element.val().length === 0) {
                    allFilled = false;
                }
            }
            return allFilled;
        };
        /**
         * Default function for getting select conditions
         */
        Criteria.inputValueSelect = function (el) {
            var values = [];
            // Go through the select elements and push each selected option to the return array
            for (var _i = 0, el_3 = el; _i < el_3.length; _i++) {
                var element = el_3[_i];
                if (element.is('select')) {
                    values.push(Criteria._escapeHTML(element.children('option:selected').data('sbv')));
                }
            }
            return values;
        };
        /**
         * Default function for getting input conditions
         */
        Criteria.inputValueInput = function (el) {
            var values = [];
            // Go through the input elements and push each value to the return array
            for (var _i = 0, el_4 = el; _i < el_4.length; _i++) {
                var element = el_4[_i];
                if (element.is('input')) {
                    values.push(Criteria._escapeHTML(element.val()));
                }
            }
            return values;
        };
        /**
         * Function that is run on each element as a call back when a search should be triggered
         */
        Criteria.updateListener = function (that, el, code) {
            // When the value is changed the criteria is now complete so can be included in searches
            // Get the condition from the map based on the key that has been selected for the condition
            var condition = that.s.conditions[that.s.condition];
            that.s.filled = condition.isInputValid(that.dom.value, that);
            that.s.value = condition.inputValue(that.dom.value, that);
            if (!that.s.filled) {
                if (!that.c.enterSearch &&
                    !(that.s.dt.settings()[0].oInit.search !== undefined &&
                        that.s.dt.settings()[0].oInit.search["return"]) ||
                    code === 13) {
                    that.s.dt.draw();
                }
                return;
            }
            if (!Array.isArray(that.s.value)) {
                that.s.value = [that.s.value];
            }
            for (var i = 0; i < that.s.value.length; i++) {
                // If the value is an array we need to sort it
                if (Array.isArray(that.s.value[i])) {
                    that.s.value[i].sort();
                }
                // Otherwise replace the decimal place character for i18n
                else if (that.s.type.includes('num') &&
                    (that.s.dt.settings()[0].oLanguage.sDecimal !== '' ||
                        that.s.dt.settings()[0].oLanguage.sThousands !== '')) {
                    var splitRD = [that.s.value[i].toString()];
                    if (that.s.dt.settings()[0].oLanguage.sDecimal !== '') {
                        splitRD = that.s.value[i].split(that.s.dt.settings()[0].oLanguage.sDecimal);
                    }
                    if (that.s.dt.settings()[0].oLanguage.sThousands !== '') {
                        for (var j = 0; j < splitRD.length; j++) {
                            splitRD[j] = splitRD[j].replace(that.s.dt.settings()[0].oLanguage.sThousands, ',');
                        }
                    }
                    that.s.value[i] = splitRD.join('.');
                }
            }
            // Take note of the cursor position so that we can refocus there later
            var idx = null;
            var cursorPos = null;
            for (var i = 0; i < that.dom.value.length; i++) {
                if (el === that.dom.value[i][0]) {
                    idx = i;
                    if (el.selectionStart !== undefined) {
                        cursorPos = el.selectionStart;
                    }
                }
            }
            if (!that.c.enterSearch &&
                !(that.s.dt.settings()[0].oInit.search !== undefined &&
                    that.s.dt.settings()[0].oInit.search["return"]) ||
                code === 13) {
                // Trigger a search
                that.s.dt.draw();
            }
            // Refocus the element and set the correct cursor position
            if (idx !== null) {
                that.dom.value[idx].removeClass(that.classes.italic);
                that.dom.value[idx].focus();
                if (cursorPos !== null) {
                    that.dom.value[idx][0].setSelectionRange(cursorPos, cursorPos);
                }
            }
        };
        // The order of the conditions will make eslint sad :(
        // Has to be in this order so that they are displayed correctly in select elements
        // Also have to disable member ordering for this as the private methods used are not yet declared otherwise
        // eslint-disable-next-line @typescript-eslint/member-ordering
        Criteria.dateConditions = {
            '=': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.date.equals', i18n.conditions.date.equals);
                },
                init: Criteria.initDate,
                inputValue: Criteria.inputValueInput,
                isInputValid: Criteria.isInputValidInput,
                search: function (value, comparison) {
                    value = value.replace(/(\/|-|,)/g, '-');
                    return value === comparison[0];
                }
            },
            // eslint-disable-next-line sort-keys
            '!=': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.date.not', i18n.conditions.date.not);
                },
                init: Criteria.initDate,
                inputValue: Criteria.inputValueInput,
                isInputValid: Criteria.isInputValidInput,
                search: function (value, comparison) {
                    value = value.replace(/(\/|-|,)/g, '-');
                    return value !== comparison[0];
                }
            },
            '<': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.date.before', i18n.conditions.date.before);
                },
                init: Criteria.initDate,
                inputValue: Criteria.inputValueInput,
                isInputValid: Criteria.isInputValidInput,
                search: function (value, comparison) {
                    value = value.replace(/(\/|-|,)/g, '-');
                    return value < comparison[0];
                }
            },
            '>': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.date.after', i18n.conditions.date.after);
                },
                init: Criteria.initDate,
                inputValue: Criteria.inputValueInput,
                isInputValid: Criteria.isInputValidInput,
                search: function (value, comparison) {
                    value = value.replace(/(\/|-|,)/g, '-');
                    return value > comparison[0];
                }
            },
            'between': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.date.between', i18n.conditions.date.between);
                },
                init: Criteria.init2Date,
                inputValue: Criteria.inputValueInput,
                isInputValid: Criteria.isInputValidInput,
                search: function (value, comparison) {
                    value = value.replace(/(\/|-|,)/g, '-');
                    if (comparison[0] < comparison[1]) {
                        return comparison[0] <= value && value <= comparison[1];
                    }
                    else {
                        return comparison[1] <= value && value <= comparison[0];
                    }
                }
            },
            // eslint-disable-next-line sort-keys
            '!between': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.date.notBetween', i18n.conditions.date.notBetween);
                },
                init: Criteria.init2Date,
                inputValue: Criteria.inputValueInput,
                isInputValid: Criteria.isInputValidInput,
                search: function (value, comparison) {
                    value = value.replace(/(\/|-|,)/g, '-');
                    if (comparison[0] < comparison[1]) {
                        return !(comparison[0] <= value && value <= comparison[1]);
                    }
                    else {
                        return !(comparison[1] <= value && value <= comparison[0]);
                    }
                }
            },
            'null': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.date.empty', i18n.conditions.date.empty);
                },
                init: Criteria.initNoValue,
                inputValue: function () {
                    return;
                },
                isInputValid: function () {
                    return true;
                },
                search: function (value) {
                    return value === null || value === undefined || value.length === 0;
                }
            },
            // eslint-disable-next-line sort-keys
            '!null': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.date.notEmpty', i18n.conditions.date.notEmpty);
                },
                init: Criteria.initNoValue,
                inputValue: function () {
                    return;
                },
                isInputValid: function () {
                    return true;
                },
                search: function (value) {
                    return !(value === null || value === undefined || value.length === 0);
                }
            }
        };
        // The order of the conditions will make eslint sad :(
        // Has to be in this order so that they are displayed correctly in select elements
        // Also have to disable member ordering for this as the private methods used are not yet declared otherwise
        // eslint-disable-next-line @typescript-eslint/member-ordering
        Criteria.momentDateConditions = {
            '=': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.date.equals', i18n.conditions.date.equals);
                },
                init: Criteria.initDate,
                inputValue: Criteria.inputValueInput,
                isInputValid: Criteria.isInputValidInput,
                search: function (value, comparison, that) {
                    return moment()(value, that.s.dateFormat).valueOf() ===
                        moment()(comparison[0], that.s.dateFormat).valueOf();
                }
            },
            // eslint-disable-next-line sort-keys
            '!=': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.date.not', i18n.conditions.date.not);
                },
                init: Criteria.initDate,
                inputValue: Criteria.inputValueInput,
                isInputValid: Criteria.isInputValidInput,
                search: function (value, comparison, that) {
                    return moment()(value, that.s.dateFormat).valueOf() !==
                        moment()(comparison[0], that.s.dateFormat).valueOf();
                }
            },
            '<': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.date.before', i18n.conditions.date.before);
                },
                init: Criteria.initDate,
                inputValue: Criteria.inputValueInput,
                isInputValid: Criteria.isInputValidInput,
                search: function (value, comparison, that) {
                    return moment()(value, that.s.dateFormat).valueOf() < moment()(comparison[0], that.s.dateFormat).valueOf();
                }
            },
            '>': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.date.after', i18n.conditions.date.after);
                },
                init: Criteria.initDate,
                inputValue: Criteria.inputValueInput,
                isInputValid: Criteria.isInputValidInput,
                search: function (value, comparison, that) {
                    return moment()(value, that.s.dateFormat).valueOf() > moment()(comparison[0], that.s.dateFormat).valueOf();
                }
            },
            'between': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.date.between', i18n.conditions.date.between);
                },
                init: Criteria.init2Date,
                inputValue: Criteria.inputValueInput,
                isInputValid: Criteria.isInputValidInput,
                search: function (value, comparison, that) {
                    var val = moment()(value, that.s.dateFormat).valueOf();
                    var comp0 = moment()(comparison[0], that.s.dateFormat).valueOf();
                    var comp1 = moment()(comparison[1], that.s.dateFormat).valueOf();
                    if (comp0 < comp1) {
                        return comp0 <= val && val <= comp1;
                    }
                    else {
                        return comp1 <= val && val <= comp0;
                    }
                }
            },
            // eslint-disable-next-line sort-keys
            '!between': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.date.notBetween', i18n.conditions.date.notBetween);
                },
                init: Criteria.init2Date,
                inputValue: Criteria.inputValueInput,
                isInputValid: Criteria.isInputValidInput,
                search: function (value, comparison, that) {
                    var val = moment()(value, that.s.dateFormat).valueOf();
                    var comp0 = moment()(comparison[0], that.s.dateFormat).valueOf();
                    var comp1 = moment()(comparison[1], that.s.dateFormat).valueOf();
                    if (comp0 < comp1) {
                        return !(+comp0 <= +val && +val <= +comp1);
                    }
                    else {
                        return !(+comp1 <= +val && +val <= +comp0);
                    }
                }
            },
            'null': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.date.empty', i18n.conditions.date.empty);
                },
                init: Criteria.initNoValue,
                inputValue: function () {
                    return;
                },
                isInputValid: function () {
                    return true;
                },
                search: function (value) {
                    return value === null || value === undefined || value.length === 0;
                }
            },
            // eslint-disable-next-line sort-keys
            '!null': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.date.notEmpty', i18n.conditions.date.notEmpty);
                },
                init: Criteria.initNoValue,
                inputValue: function () {
                    return;
                },
                isInputValid: function () {
                    return true;
                },
                search: function (value) {
                    return !(value === null || value === undefined || value.length === 0);
                }
            }
        };
        // The order of the conditions will make eslint sad :(
        // Has to be in this order so that they are displayed correctly in select elements
        // Also have to disable member ordering for this as the private methods used are not yet declared otherwise
        // eslint-disable-next-line @typescript-eslint/member-ordering
        Criteria.luxonDateConditions = {
            '=': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.date.equals', i18n.conditions.date.equals);
                },
                init: Criteria.initDate,
                inputValue: Criteria.inputValueInput,
                isInputValid: Criteria.isInputValidInput,
                search: function (value, comparison, that) {
                    return luxon().DateTime.fromFormat(value, that.s.dateFormat).ts
                        === luxon().DateTime.fromFormat(comparison[0], that.s.dateFormat).ts;
                }
            },
            // eslint-disable-next-line sort-keys
            '!=': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.date.not', i18n.conditions.date.not);
                },
                init: Criteria.initDate,
                inputValue: Criteria.inputValueInput,
                isInputValid: Criteria.isInputValidInput,
                search: function (value, comparison, that) {
                    return luxon().DateTime.fromFormat(value, that.s.dateFormat).ts
                        !== luxon().DateTime.fromFormat(comparison[0], that.s.dateFormat).ts;
                }
            },
            '<': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.date.before', i18n.conditions.date.before);
                },
                init: Criteria.initDate,
                inputValue: Criteria.inputValueInput,
                isInputValid: Criteria.isInputValidInput,
                search: function (value, comparison, that) {
                    return luxon().DateTime.fromFormat(value, that.s.dateFormat).ts
                        < luxon().DateTime.fromFormat(comparison[0], that.s.dateFormat).ts;
                }
            },
            '>': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.date.after', i18n.conditions.date.after);
                },
                init: Criteria.initDate,
                inputValue: Criteria.inputValueInput,
                isInputValid: Criteria.isInputValidInput,
                search: function (value, comparison, that) {
                    return luxon().DateTime.fromFormat(value, that.s.dateFormat).ts
                        > luxon().DateTime.fromFormat(comparison[0], that.s.dateFormat).ts;
                }
            },
            'between': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.date.between', i18n.conditions.date.between);
                },
                init: Criteria.init2Date,
                inputValue: Criteria.inputValueInput,
                isInputValid: Criteria.isInputValidInput,
                search: function (value, comparison, that) {
                    var val = luxon().DateTime.fromFormat(value, that.s.dateFormat).ts;
                    var comp0 = luxon().DateTime.fromFormat(comparison[0], that.s.dateFormat).ts;
                    var comp1 = luxon().DateTime.fromFormat(comparison[1], that.s.dateFormat).ts;
                    if (comp0 < comp1) {
                        return comp0 <= val && val <= comp1;
                    }
                    else {
                        return comp1 <= val && val <= comp0;
                    }
                }
            },
            // eslint-disable-next-line sort-keys
            '!between': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.date.notBetween', i18n.conditions.date.notBetween);
                },
                init: Criteria.init2Date,
                inputValue: Criteria.inputValueInput,
                isInputValid: Criteria.isInputValidInput,
                search: function (value, comparison, that) {
                    var val = luxon().DateTime.fromFormat(value, that.s.dateFormat).ts;
                    var comp0 = luxon().DateTime.fromFormat(comparison[0], that.s.dateFormat).ts;
                    var comp1 = luxon().DateTime.fromFormat(comparison[1], that.s.dateFormat).ts;
                    if (comp0 < comp1) {
                        return !(+comp0 <= +val && +val <= +comp1);
                    }
                    else {
                        return !(+comp1 <= +val && +val <= +comp0);
                    }
                }
            },
            'null': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.date.empty', i18n.conditions.date.empty);
                },
                init: Criteria.initNoValue,
                inputValue: function () {
                    return;
                },
                isInputValid: function () {
                    return true;
                },
                search: function (value) {
                    return value === null || value === undefined || value.length === 0;
                }
            },
            // eslint-disable-next-line sort-keys
            '!null': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.date.notEmpty', i18n.conditions.date.notEmpty);
                },
                init: Criteria.initNoValue,
                inputValue: function () {
                    return;
                },
                isInputValid: function () {
                    return true;
                },
                search: function (value) {
                    return !(value === null || value === undefined || value.length === 0);
                }
            }
        };
        // The order of the conditions will make eslint sad :(
        // Has to be in this order so that they are displayed correctly in select elements
        // Also have to disable member ordering for this as the private methods used are not yet declared otherwise
        // eslint-disable-next-line @typescript-eslint/member-ordering
        Criteria.numConditions = {
            '=': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.number.equals', i18n.conditions.number.equals);
                },
                init: Criteria.initSelect,
                inputValue: Criteria.inputValueSelect,
                isInputValid: Criteria.isInputValidSelect,
                search: function (value, comparison) {
                    return +value === +comparison[0];
                }
            },
            // eslint-disable-next-line sort-keys
            '!=': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.number.not', i18n.conditions.number.not);
                },
                init: Criteria.initSelect,
                inputValue: Criteria.inputValueSelect,
                isInputValid: Criteria.isInputValidSelect,
                search: function (value, comparison) {
                    return +value !== +comparison[0];
                }
            },
            '<': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.number.lt', i18n.conditions.number.lt);
                },
                init: Criteria.initInput,
                inputValue: Criteria.inputValueInput,
                isInputValid: Criteria.isInputValidInput,
                search: function (value, comparison) {
                    return +value < +comparison[0];
                }
            },
            '<=': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.number.lte', i18n.conditions.number.lte);
                },
                init: Criteria.initInput,
                inputValue: Criteria.inputValueInput,
                isInputValid: Criteria.isInputValidInput,
                search: function (value, comparison) {
                    return +value <= +comparison[0];
                }
            },
            '>=': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.number.gte', i18n.conditions.number.gte);
                },
                init: Criteria.initInput,
                inputValue: Criteria.inputValueInput,
                isInputValid: Criteria.isInputValidInput,
                search: function (value, comparison) {
                    return +value >= +comparison[0];
                }
            },
            // eslint-disable-next-line sort-keys
            '>': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.number.gt', i18n.conditions.number.gt);
                },
                init: Criteria.initInput,
                inputValue: Criteria.inputValueInput,
                isInputValid: Criteria.isInputValidInput,
                search: function (value, comparison) {
                    return +value > +comparison[0];
                }
            },
            'between': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.number.between', i18n.conditions.number.between);
                },
                init: Criteria.init2Input,
                inputValue: Criteria.inputValueInput,
                isInputValid: Criteria.isInputValidInput,
                search: function (value, comparison) {
                    if (+comparison[0] < +comparison[1]) {
                        return +comparison[0] <= +value && +value <= +comparison[1];
                    }
                    else {
                        return +comparison[1] <= +value && +value <= +comparison[0];
                    }
                }
            },
            // eslint-disable-next-line sort-keys
            '!between': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.number.notBetween', i18n.conditions.number.notBetween);
                },
                init: Criteria.init2Input,
                inputValue: Criteria.inputValueInput,
                isInputValid: Criteria.isInputValidInput,
                search: function (value, comparison) {
                    if (+comparison[0] < +comparison[1]) {
                        return !(+comparison[0] <= +value && +value <= +comparison[1]);
                    }
                    else {
                        return !(+comparison[1] <= +value && +value <= +comparison[0]);
                    }
                }
            },
            'null': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.number.empty', i18n.conditions.number.empty);
                },
                init: Criteria.initNoValue,
                inputValue: function () {
                    return;
                },
                isInputValid: function () {
                    return true;
                },
                search: function (value) {
                    return value === null || value === undefined || value.length === 0;
                }
            },
            // eslint-disable-next-line sort-keys
            '!null': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.number.notEmpty', i18n.conditions.number.notEmpty);
                },
                init: Criteria.initNoValue,
                inputValue: function () {
                    return;
                },
                isInputValid: function () {
                    return true;
                },
                search: function (value) {
                    return !(value === null || value === undefined || value.length === 0);
                }
            }
        };
        // The order of the conditions will make eslint sad :(
        // Has to be in this order so that they are displayed correctly in select elements
        // Also have to disable member ordering for this as the private methods used are not yet declared otherwise
        // eslint-disable-next-line @typescript-eslint/member-ordering
        Criteria.numFmtConditions = {
            '=': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.number.equals', i18n.conditions.number.equals);
                },
                init: Criteria.initSelect,
                inputValue: Criteria.inputValueSelect,
                isInputValid: Criteria.isInputValidSelect,
                search: function (value, comparison) {
                    return Criteria.parseNumFmt(value) === Criteria.parseNumFmt(comparison[0]);
                }
            },
            // eslint-disable-next-line sort-keys
            '!=': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.number.not', i18n.conditions.number.not);
                },
                init: Criteria.initSelect,
                inputValue: Criteria.inputValueSelect,
                isInputValid: Criteria.isInputValidSelect,
                search: function (value, comparison) {
                    return Criteria.parseNumFmt(value) !== Criteria.parseNumFmt(comparison[0]);
                }
            },
            '<': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.number.lt', i18n.conditions.number.lt);
                },
                init: Criteria.initInput,
                inputValue: Criteria.inputValueInput,
                isInputValid: Criteria.isInputValidInput,
                search: function (value, comparison) {
                    return Criteria.parseNumFmt(value) < Criteria.parseNumFmt(comparison[0]);
                }
            },
            '<=': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.number.lte', i18n.conditions.number.lte);
                },
                init: Criteria.initInput,
                inputValue: Criteria.inputValueInput,
                isInputValid: Criteria.isInputValidInput,
                search: function (value, comparison) {
                    return Criteria.parseNumFmt(value) <= Criteria.parseNumFmt(comparison[0]);
                }
            },
            '>=': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.number.gte', i18n.conditions.number.gte);
                },
                init: Criteria.initInput,
                inputValue: Criteria.inputValueInput,
                isInputValid: Criteria.isInputValidInput,
                search: function (value, comparison) {
                    return Criteria.parseNumFmt(value) >= Criteria.parseNumFmt(comparison[0]);
                }
            },
            // eslint-disable-next-line sort-keys
            '>': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.number.gt', i18n.conditions.number.gt);
                },
                init: Criteria.initInput,
                inputValue: Criteria.inputValueInput,
                isInputValid: Criteria.isInputValidInput,
                search: function (value, comparison) {
                    return Criteria.parseNumFmt(value) > Criteria.parseNumFmt(comparison[0]);
                }
            },
            'between': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.number.between', i18n.conditions.number.between);
                },
                init: Criteria.init2Input,
                inputValue: Criteria.inputValueInput,
                isInputValid: Criteria.isInputValidInput,
                search: function (value, comparison) {
                    var val = Criteria.parseNumFmt(value);
                    var comp0 = Criteria.parseNumFmt(comparison[0]);
                    var comp1 = Criteria.parseNumFmt(comparison[1]);
                    if (+comp0 < +comp1) {
                        return +comp0 <= +val && +val <= +comp1;
                    }
                    else {
                        return +comp1 <= +val && +val <= +comp0;
                    }
                }
            },
            // eslint-disable-next-line sort-keys
            '!between': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.number.notBetween', i18n.conditions.number.notBetween);
                },
                init: Criteria.init2Input,
                inputValue: Criteria.inputValueInput,
                isInputValid: Criteria.isInputValidInput,
                search: function (value, comparison) {
                    var val = Criteria.parseNumFmt(value);
                    var comp0 = Criteria.parseNumFmt(comparison[0]);
                    var comp1 = Criteria.parseNumFmt(comparison[1]);
                    if (+comp0 < +comp1) {
                        return !(+comp0 <= +val && +val <= +comp1);
                    }
                    else {
                        return !(+comp1 <= +val && +val <= +comp0);
                    }
                }
            },
            'null': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.number.empty', i18n.conditions.number.empty);
                },
                init: Criteria.initNoValue,
                inputValue: function () {
                    return;
                },
                isInputValid: function () {
                    return true;
                },
                search: function (value) {
                    return value === null || value === undefined || value.length === 0;
                }
            },
            // eslint-disable-next-line sort-keys
            '!null': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.number.notEmpty', i18n.conditions.number.notEmpty);
                },
                init: Criteria.initNoValue,
                inputValue: function () {
                    return;
                },
                isInputValid: function () {
                    return true;
                },
                search: function (value) {
                    return !(value === null || value === undefined || value.length === 0);
                }
            }
        };
        // The order of the conditions will make eslint sad :(
        // Has to be in this order so that they are displayed correctly in select elements
        // Also have to disable member ordering for this as the private methods used are not yet declared otherwise
        // eslint-disable-next-line @typescript-eslint/member-ordering
        Criteria.stringConditions = {
            '=': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.string.equals', i18n.conditions.string.equals);
                },
                init: Criteria.initSelect,
                inputValue: Criteria.inputValueSelect,
                isInputValid: Criteria.isInputValidSelect,
                search: function (value, comparison) {
                    return value === comparison[0];
                }
            },
            // eslint-disable-next-line sort-keys
            '!=': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.string.not', i18n.conditions.string.not);
                },
                init: Criteria.initSelect,
                inputValue: Criteria.inputValueSelect,
                isInputValid: Criteria.isInputValidInput,
                search: function (value, comparison) {
                    return value !== comparison[0];
                }
            },
            'starts': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.string.startsWith', i18n.conditions.string.startsWith);
                },
                init: Criteria.initInput,
                inputValue: Criteria.inputValueInput,
                isInputValid: Criteria.isInputValidInput,
                search: function (value, comparison) {
                    return value.toLowerCase().indexOf(comparison[0].toLowerCase()) === 0;
                }
            },
            // eslint-disable-next-line sort-keys
            '!starts': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.string.notStartsWith', i18n.conditions.string.notStartsWith);
                },
                init: Criteria.initInput,
                inputValue: Criteria.inputValueInput,
                isInputValid: Criteria.isInputValidInput,
                search: function (value, comparison) {
                    return value.toLowerCase().indexOf(comparison[0].toLowerCase()) !== 0;
                }
            },
            // eslint-disable-next-line sort-keys
            'contains': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.string.contains', i18n.conditions.string.contains);
                },
                init: Criteria.initInput,
                inputValue: Criteria.inputValueInput,
                isInputValid: Criteria.isInputValidInput,
                search: function (value, comparison) {
                    return value.toLowerCase().includes(comparison[0].toLowerCase());
                }
            },
            // eslint-disable-next-line sort-keys
            '!contains': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.string.notContains', i18n.conditions.string.notContains);
                },
                init: Criteria.initInput,
                inputValue: Criteria.inputValueInput,
                isInputValid: Criteria.isInputValidInput,
                search: function (value, comparison) {
                    return !value.toLowerCase().includes(comparison[0].toLowerCase());
                }
            },
            'ends': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.string.endsWith', i18n.conditions.string.endsWith);
                },
                init: Criteria.initInput,
                inputValue: Criteria.inputValueInput,
                isInputValid: Criteria.isInputValidInput,
                search: function (value, comparison) {
                    return value.toLowerCase().endsWith(comparison[0].toLowerCase());
                }
            },
            // eslint-disable-next-line sort-keys
            '!ends': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.string.notEndsWith', i18n.conditions.string.notEndsWith);
                },
                init: Criteria.initInput,
                inputValue: Criteria.inputValueInput,
                isInputValid: Criteria.isInputValidInput,
                search: function (value, comparison) {
                    return !value.toLowerCase().endsWith(comparison[0].toLowerCase());
                }
            },
            'null': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.string.empty', i18n.conditions.string.empty);
                },
                init: Criteria.initNoValue,
                inputValue: function () {
                    return;
                },
                isInputValid: function () {
                    return true;
                },
                search: function (value) {
                    return value === null || value === undefined || value.length === 0;
                }
            },
            // eslint-disable-next-line sort-keys
            '!null': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.string.notEmpty', i18n.conditions.string.notEmpty);
                },
                init: Criteria.initNoValue,
                inputValue: function () {
                    return;
                },
                isInputValid: function () {
                    return true;
                },
                search: function (value) {
                    return !(value === null || value === undefined || value.length === 0);
                }
            }
        };
        // The order of the conditions will make eslint sad :(
        // Also have to disable member ordering for this as the private methods used are not yet declared otherwise
        // eslint-disable-next-line @typescript-eslint/member-ordering
        Criteria.arrayConditions = {
            'contains': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.array.contains', i18n.conditions.array.contains);
                },
                init: Criteria.initSelectArray,
                inputValue: Criteria.inputValueSelect,
                isInputValid: Criteria.isInputValidSelect,
                search: function (value, comparison) {
                    return value.includes(comparison[0]);
                }
            },
            'without': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.array.without', i18n.conditions.array.without);
                },
                init: Criteria.initSelectArray,
                inputValue: Criteria.inputValueSelect,
                isInputValid: Criteria.isInputValidSelect,
                search: function (value, comparison) {
                    return value.indexOf(comparison[0]) === -1;
                }
            },
            // eslint-disable-next-line sort-keys
            '=': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.array.equals', i18n.conditions.array.equals);
                },
                init: Criteria.initSelect,
                inputValue: Criteria.inputValueSelect,
                isInputValid: Criteria.isInputValidSelect,
                search: function (value, comparison) {
                    if (value.length === comparison[0].length) {
                        for (var i = 0; i < value.length; i++) {
                            if (value[i] !== comparison[0][i]) {
                                return false;
                            }
                        }
                        return true;
                    }
                    return false;
                }
            },
            // eslint-disable-next-line sort-keys
            '!=': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.array.not', i18n.conditions.array.not);
                },
                init: Criteria.initSelect,
                inputValue: Criteria.inputValueSelect,
                isInputValid: Criteria.isInputValidSelect,
                search: function (value, comparison) {
                    if (value.length === comparison[0].length) {
                        for (var i = 0; i < value.length; i++) {
                            if (value[i] !== comparison[0][i]) {
                                return true;
                            }
                        }
                        return false;
                    }
                    return true;
                }
            },
            'null': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.array.empty', i18n.conditions.array.empty);
                },
                init: Criteria.initNoValue,
                inputValue: function () {
                    return;
                },
                isInputValid: function () {
                    return true;
                },
                search: function (value) {
                    return value === null || value === undefined || value.length === 0;
                }
            },
            // eslint-disable-next-line sort-keys
            '!null': {
                conditionName: function (dt, i18n) {
                    return dt.i18n('searchBuilder.conditions.array.notEmpty', i18n.conditions.array.notEmpty);
                },
                init: Criteria.initNoValue,
                inputValue: function () {
                    return;
                },
                isInputValid: function () {
                    return true;
                },
                search: function (value) {
                    return value !== null && value !== undefined && value.length !== 0;
                }
            }
        };
        // eslint will be sad because we have to disable member ordering for this as the
        // private static properties used are not yet declared otherwise
        // eslint-disable-next-line @typescript-eslint/member-ordering
        Criteria.defaults = {
            columns: true,
            conditions: {
                'array': Criteria.arrayConditions,
                'date': Criteria.dateConditions,
                'html': Criteria.stringConditions,
                'html-num': Criteria.numConditions,
                'html-num-fmt': Criteria.numFmtConditions,
                'luxon': Criteria.luxonDateConditions,
                'moment': Criteria.momentDateConditions,
                'num': Criteria.numConditions,
                'num-fmt': Criteria.numFmtConditions,
                'string': Criteria.stringConditions
            },
            depthLimit: false,
            enterSearch: false,
            filterChanged: undefined,
            greyscale: false,
            i18n: {
                add: 'Add Condition',
                button: {
                    0: 'Search Builder',
                    _: 'Search Builder (%d)'
                },
                clearAll: 'Clear All',
                condition: 'Condition',
                data: 'Data',
                "delete": '&times',
                deleteTitle: 'Delete filtering rule',
                left: '<',
                leftTitle: 'Outdent criteria',
                logicAnd: 'And',
                logicOr: 'Or',
                right: '>',
                rightTitle: 'Indent criteria',
                title: {
                    0: 'Custom Search Builder',
                    _: 'Custom Search Builder (%d)'
                },
                value: 'Value',
                valueJoiner: 'and'
            },
            logic: 'AND',
            orthogonal: {
                display: 'display',
                search: 'filter'
            },
            preDefined: false
        };
        return Criteria;
    }());

    var $$2;
    var dataTable$2;
    /**
     * Sets the value of jQuery for use in the file
     *
     * @param jq the instance of jQuery to be set
     */
    function setJQuery$1(jq) {
        $$2 = jq;
        dataTable$2 = jq.fn.dataTable;
    }
    /**
     * The Group class is used within SearchBuilder to represent a group of criteria
     */
    var Group = /** @class */ (function () {
        function Group(table, opts, topGroup, index, isChild, depth, serverData) {
            if (index === void 0) { index = 0; }
            if (isChild === void 0) { isChild = false; }
            if (depth === void 0) { depth = 1; }
            if (serverData === void 0) { serverData = undefined; }
            // Check that the required version of DataTables is included
            if (!dataTable$2 || !dataTable$2.versionCheck || !dataTable$2.versionCheck('1.10.0')) {
                throw new Error('SearchBuilder requires DataTables 1.10 or newer');
            }
            this.classes = $$2.extend(true, {}, Group.classes);
            // Get options from user
            this.c = $$2.extend(true, {}, Group.defaults, opts);
            this.s = {
                criteria: [],
                depth: depth,
                dt: table,
                index: index,
                isChild: isChild,
                logic: undefined,
                opts: opts,
                preventRedraw: false,
                serverData: serverData,
                toDrop: undefined,
                topGroup: topGroup
            };
            this.dom = {
                add: $$2('<button/>')
                    .addClass(this.classes.add)
                    .addClass(this.classes.button)
                    .attr('type', 'button'),
                clear: $$2('<button>&times</button>')
                    .addClass(this.classes.button)
                    .addClass(this.classes.clearGroup)
                    .attr('type', 'button'),
                container: $$2('<div/>')
                    .addClass(this.classes.group),
                logic: $$2('<button><div/></button>')
                    .addClass(this.classes.logic)
                    .addClass(this.classes.button)
                    .attr('type', 'button'),
                logicContainer: $$2('<div/>')
                    .addClass(this.classes.logicContainer)
            };
            // A reference to the top level group is maintained throughout any subgroups and criteria that may be created
            if (this.s.topGroup === undefined) {
                this.s.topGroup = this.dom.container;
            }
            this._setup();
            return this;
        }
        /**
         * Destroys the groups buttons, clears the internal criteria and removes it from the dom
         */
        Group.prototype.destroy = function () {
            // Turn off listeners
            this.dom.add.off('.dtsb');
            this.dom.logic.off('.dtsb');
            // Trigger event for groups at a higher level to pick up on
            this.dom.container
                .trigger('dtsb-destroy')
                .remove();
            this.s.criteria = [];
        };
        /**
         * Gets the details required to rebuild the group
         */
        // Eslint upset at empty object but needs to be done
        // eslint-disable-next-line @typescript-eslint/ban-types
        Group.prototype.getDetails = function (deFormatDates) {
            if (deFormatDates === void 0) { deFormatDates = false; }
            if (this.s.criteria.length === 0) {
                return {};
            }
            var details = {
                criteria: [],
                logic: this.s.logic
            };
            // NOTE here crit could be either a subgroup or a criteria
            for (var _i = 0, _a = this.s.criteria; _i < _a.length; _i++) {
                var crit = _a[_i];
                details.criteria.push(crit.criteria.getDetails(deFormatDates));
            }
            return details;
        };
        /**
         * Getter for the node for the container of the group
         *
         * @returns Node for the container of the group
         */
        Group.prototype.getNode = function () {
            return this.dom.container;
        };
        /**
         * Rebuilds the group based upon the details passed in
         *
         * @param loadedDetails the details required to rebuild the group
         */
        Group.prototype.rebuild = function (loadedDetails) {
            // If no criteria are stored then just return
            if (loadedDetails.criteria === undefined ||
                loadedDetails.criteria === null ||
                Array.isArray(loadedDetails.criteria) && loadedDetails.criteria.length === 0) {
                return;
            }
            this.s.logic = loadedDetails.logic;
            this.dom.logic.children().first().html(this.s.logic === 'OR'
                ? this.s.dt.i18n('searchBuilder.logicOr', this.c.i18n.logicOr)
                : this.s.dt.i18n('searchBuilder.logicAnd', this.c.i18n.logicAnd));
            // Add all of the criteria, be it a sub group or a criteria
            if (Array.isArray(loadedDetails.criteria)) {
                for (var _i = 0, _a = loadedDetails.criteria; _i < _a.length; _i++) {
                    var crit = _a[_i];
                    if (crit.logic !== undefined) {
                        this._addPrevGroup(crit);
                    }
                    else if (crit.logic === undefined) {
                        this._addPrevCriteria(crit);
                    }
                }
            }
            // For all of the criteria children, update the arrows incase they require changing and set the listeners
            for (var _b = 0, _c = this.s.criteria; _b < _c.length; _b++) {
                var crit = _c[_b];
                if (crit.criteria instanceof Criteria) {
                    crit.criteria.updateArrows(this.s.criteria.length > 1);
                    this._setCriteriaListeners(crit.criteria);
                }
            }
        };
        /**
         * Redraws the Contents of the searchBuilder Groups and Criteria
         */
        Group.prototype.redrawContents = function () {
            if (this.s.preventRedraw) {
                return;
            }
            // Clear the container out and add the basic elements
            this.dom.container.children().detach();
            this.dom.container
                .append(this.dom.logicContainer)
                .append(this.dom.add);
            // Sort the criteria by index so that they appear in the correct order
            this.s.criteria.sort(function (a, b) {
                if (a.criteria.s.index < b.criteria.s.index) {
                    return -1;
                }
                else if (a.criteria.s.index > b.criteria.s.index) {
                    return 1;
                }
                return 0;
            });
            this.setListeners();
            for (var i = 0; i < this.s.criteria.length; i++) {
                var crit = this.s.criteria[i].criteria;
                if (crit instanceof Criteria) {
                    // Reset the index to the new value
                    this.s.criteria[i].index = i;
                    this.s.criteria[i].criteria.s.index = i;
                    // Add to the group
                    this.s.criteria[i].criteria.dom.container.insertBefore(this.dom.add);
                    // Set listeners for various points
                    this._setCriteriaListeners(crit);
                    this.s.criteria[i].criteria.s.preventRedraw = this.s.preventRedraw;
                    this.s.criteria[i].criteria.rebuild(this.s.criteria[i].criteria.getDetails());
                    this.s.criteria[i].criteria.s.preventRedraw = false;
                }
                else if (crit instanceof Group && crit.s.criteria.length > 0) {
                    // Reset the index to the new value
                    this.s.criteria[i].index = i;
                    this.s.criteria[i].criteria.s.index = i;
                    // Add the sub group to the group
                    this.s.criteria[i].criteria.dom.container.insertBefore(this.dom.add);
                    // Redraw the contents of the group
                    crit.s.preventRedraw = this.s.preventRedraw;
                    crit.redrawContents();
                    crit.s.preventRedraw = false;
                    this._setGroupListeners(crit);
                }
                else {
                    // The group is empty so remove it
                    this.s.criteria.splice(i, 1);
                    i--;
                }
            }
            this.setupLogic();
        };
        /**
         * Resizes the logic button only rather than the entire dom.
         */
        Group.prototype.redrawLogic = function () {
            for (var _i = 0, _a = this.s.criteria; _i < _a.length; _i++) {
                var crit = _a[_i];
                if (crit.criteria instanceof Group) {
                    crit.criteria.redrawLogic();
                }
            }
            this.setupLogic();
        };
        /**
         * Search method, checking the row data against the criteria in the group
         *
         * @param rowData The row data to be compared
         * @returns boolean The result of the search
         */
        Group.prototype.search = function (rowData, rowIdx) {
            if (this.s.logic === 'AND') {
                return this._andSearch(rowData, rowIdx);
            }
            else if (this.s.logic === 'OR') {
                return this._orSearch(rowData, rowIdx);
            }
            return true;
        };
        /**
         * Locates the groups logic button to the correct location on the page
         */
        Group.prototype.setupLogic = function () {
            // Remove logic button
            this.dom.logicContainer.remove();
            this.dom.clear.remove();
            // If there are no criteria in the group then keep the logic removed and return
            if (this.s.criteria.length < 1) {
                if (!this.s.isChild) {
                    this.dom.container.trigger('dtsb-destroy');
                    // Set criteria left margin
                    this.dom.container.css('margin-left', 0);
                }
                return;
            }
            this.dom.clear.height('0px');
            this.dom.logicContainer.append(this.dom.clear);
            // Prepend logic button
            this.dom.container.prepend(this.dom.logicContainer);
            for (var _i = 0, _a = this.s.criteria; _i < _a.length; _i++) {
                var crit = _a[_i];
                if (crit.criteria instanceof Criteria) {
                    crit.criteria.setupButtons();
                }
            }
            // Set width, take 2 for the border
            var height = this.dom.container.outerHeight() - 1;
            this.dom.logicContainer.width(height);
            this._setLogicListener();
            // Set criteria left margin
            this.dom.container.css('margin-left', this.dom.logicContainer.outerHeight(true));
            var logicOffset = this.dom.logicContainer.offset();
            // Set horizontal alignment
            var currentLeft = logicOffset.left;
            var groupLeft = this.dom.container.offset().left;
            var shuffleLeft = currentLeft - groupLeft;
            var newPos = currentLeft - shuffleLeft - this.dom.logicContainer.outerHeight(true);
            this.dom.logicContainer.offset({ left: newPos });
            // Set vertical alignment
            var firstCrit = this.dom.logicContainer.next();
            var currentTop = logicOffset.top;
            var firstTop = $$2(firstCrit).offset().top;
            var shuffleTop = currentTop - firstTop;
            var newTop = currentTop - shuffleTop;
            this.dom.logicContainer.offset({ top: newTop });
            this.dom.clear.outerHeight(this.dom.logicContainer.height());
            this._setClearListener();
        };
        /**
         * Sets listeners on the groups elements
         */
        Group.prototype.setListeners = function () {
            var _this = this;
            this.dom.add.unbind('click');
            this.dom.add.on('click.dtsb', function () {
                // If this is the parent group then the logic button has not been added yet
                if (!_this.s.isChild) {
                    _this.dom.container.prepend(_this.dom.logicContainer);
                }
                _this.addCriteria();
                _this.dom.container.trigger('dtsb-add');
                _this.s.dt.state.save();
                return false;
            });
            for (var _i = 0, _a = this.s.criteria; _i < _a.length; _i++) {
                var crit = _a[_i];
                crit.criteria.setListeners();
            }
            this._setClearListener();
            this._setLogicListener();
        };
        /**
         * Adds a criteria to the group
         *
         * @param crit Instance of Criteria to be added to the group
         */
        Group.prototype.addCriteria = function (crit) {
            if (crit === void 0) { crit = null; }
            var index = crit === null ? this.s.criteria.length : crit.s.index;
            var criteria = new Criteria(this.s.dt, this.s.opts, this.s.topGroup, index, this.s.depth, this.s.serverData);
            // If a Criteria has been passed in then set the values to continue that
            if (crit !== null) {
                criteria.c = crit.c;
                criteria.s = crit.s;
                criteria.s.depth = this.s.depth;
                criteria.classes = crit.classes;
            }
            criteria.populate();
            var inserted = false;
            for (var i = 0; i < this.s.criteria.length; i++) {
                if (i === 0 && this.s.criteria[i].criteria.s.index > criteria.s.index) {
                    // Add the node for the criteria at the start of the group
                    criteria.getNode().insertBefore(this.s.criteria[i].criteria.dom.container);
                    inserted = true;
                }
                else if (i < this.s.criteria.length - 1 &&
                    this.s.criteria[i].criteria.s.index < criteria.s.index &&
                    this.s.criteria[i + 1].criteria.s.index > criteria.s.index) {
                    // Add the node for the criteria in the correct location
                    criteria.getNode().insertAfter(this.s.criteria[i].criteria.dom.container);
                    inserted = true;
                }
            }
            if (!inserted) {
                criteria.getNode().insertBefore(this.dom.add);
            }
            // Add the details for this criteria to the array
            this.s.criteria.push({
                criteria: criteria,
                index: index
            });
            this.s.criteria = this.s.criteria.sort(function (a, b) { return a.criteria.s.index - b.criteria.s.index; });
            for (var _i = 0, _a = this.s.criteria; _i < _a.length; _i++) {
                var opt = _a[_i];
                if (opt.criteria instanceof Criteria) {
                    opt.criteria.updateArrows(this.s.criteria.length > 1);
                }
            }
            this._setCriteriaListeners(criteria);
            criteria.setListeners();
            this.setupLogic();
        };
        /**
         * Checks the group to see if it has any filled criteria
         */
        Group.prototype.checkFilled = function () {
            for (var _i = 0, _a = this.s.criteria; _i < _a.length; _i++) {
                var crit = _a[_i];
                if (crit.criteria instanceof Criteria && crit.criteria.s.filled ||
                    crit.criteria instanceof Group && crit.criteria.checkFilled()) {
                    return true;
                }
            }
            return false;
        };
        /**
         * Gets the count for the number of criteria in this group and any sub groups
         */
        Group.prototype.count = function () {
            var count = 0;
            for (var _i = 0, _a = this.s.criteria; _i < _a.length; _i++) {
                var crit = _a[_i];
                if (crit.criteria instanceof Group) {
                    count += crit.criteria.count();
                }
                else {
                    count++;
                }
            }
            return count;
        };
        /**
         * Rebuilds a sub group that previously existed
         *
         * @param loadedGroup The details of a group within this group
         */
        Group.prototype._addPrevGroup = function (loadedGroup) {
            var idx = this.s.criteria.length;
            var group = new Group(this.s.dt, this.c, this.s.topGroup, idx, true, this.s.depth + 1, this.s.serverData);
            // Add the new group to the criteria array
            this.s.criteria.push({
                criteria: group,
                index: idx,
                logic: group.s.logic
            });
            // Rebuild it with the previous conditions for that group
            group.rebuild(loadedGroup);
            this.s.criteria[idx].criteria = group;
            this.s.topGroup.trigger('dtsb-redrawContents');
            this._setGroupListeners(group);
        };
        /**
         * Rebuilds a criteria of this group that previously existed
         *
         * @param loadedCriteria The details of a criteria within the group
         */
        Group.prototype._addPrevCriteria = function (loadedCriteria) {
            var idx = this.s.criteria.length;
            var criteria = new Criteria(this.s.dt, this.s.opts, this.s.topGroup, idx, this.s.depth, this.s.serverData);
            criteria.populate();
            // Add the new criteria to the criteria array
            this.s.criteria.push({
                criteria: criteria,
                index: idx
            });
            // Rebuild it with the previous conditions for that criteria
            criteria.s.preventRedraw = this.s.preventRedraw;
            criteria.rebuild(loadedCriteria);
            criteria.s.preventRedraw = false;
            this.s.criteria[idx].criteria = criteria;
            if (!this.s.preventRedraw) {
                this.s.topGroup.trigger('dtsb-redrawContents');
            }
        };
        /**
         * Checks And the criteria using AND logic
         *
         * @param rowData The row data to be checked against the search criteria
         * @returns boolean The result of the AND search
         */
        Group.prototype._andSearch = function (rowData, rowIdx) {
            // If there are no criteria then return true for this group
            if (this.s.criteria.length === 0) {
                return true;
            }
            for (var _i = 0, _a = this.s.criteria; _i < _a.length; _i++) {
                var crit = _a[_i];
                // If the criteria is not complete then skip it
                if (crit.criteria instanceof Criteria && !crit.criteria.s.filled) {
                    continue;
                }
                // Otherwise if a single one fails return false
                else if (!crit.criteria.search(rowData, rowIdx)) {
                    return false;
                }
            }
            // If we get to here then everything has passed, so return true for the group
            return true;
        };
        /**
         * Checks And the criteria using OR logic
         *
         * @param rowData The row data to be checked against the search criteria
         * @returns boolean The result of the OR search
         */
        Group.prototype._orSearch = function (rowData, rowIdx) {
            // If there are no criteria in the group then return true
            if (this.s.criteria.length === 0) {
                return true;
            }
            // This will check to make sure that at least one criteria in the group is complete
            var filledfound = false;
            for (var _i = 0, _a = this.s.criteria; _i < _a.length; _i++) {
                var crit = _a[_i];
                if (crit.criteria instanceof Criteria && crit.criteria.s.filled) {
                    // A completed criteria has been found so set the flag
                    filledfound = true;
                    // If the search passes then return true
                    if (crit.criteria.search(rowData, rowIdx)) {
                        return true;
                    }
                }
                else if (crit.criteria instanceof Group && crit.criteria.checkFilled()) {
                    filledfound = true;
                    if (crit.criteria.search(rowData, rowIdx)) {
                        return true;
                    }
                }
            }
            // If we get here we need to return the inverse of filledfound,
            //  as if any have been found and we are here then none have passed
            return !filledfound;
        };
        /**
         * Removes a criteria from the group
         *
         * @param criteria The criteria instance to be removed
         */
        Group.prototype._removeCriteria = function (criteria, group) {
            if (group === void 0) { group = false; }
            // If removing a criteria and there is only then then just destroy the group
            if (this.s.criteria.length <= 1 && this.s.isChild) {
                this.destroy();
            }
            else {
                // Otherwise splice the given criteria out and redo the indexes
                var last = void 0;
                for (var i = 0; i < this.s.criteria.length; i++) {
                    if (this.s.criteria[i].index === criteria.s.index &&
                        (!group || this.s.criteria[i].criteria instanceof Group)) {
                        last = i;
                    }
                }
                // We want to remove the last element with the desired index, as its replacement will be inserted before it
                if (last !== undefined) {
                    this.s.criteria.splice(last, 1);
                }
                for (var i = 0; i < this.s.criteria.length; i++) {
                    this.s.criteria[i].index = i;
                    this.s.criteria[i].criteria.s.index = i;
                }
            }
        };
        /**
         * Sets the listeners in group for a criteria
         *
         * @param criteria The criteria for the listeners to be set on
         */
        Group.prototype._setCriteriaListeners = function (criteria) {
            var _this = this;
            criteria.dom["delete"]
                .unbind('click')
                .on('click.dtsb', function () {
                _this._removeCriteria(criteria);
                criteria.dom.container.remove();
                for (var _i = 0, _a = _this.s.criteria; _i < _a.length; _i++) {
                    var crit = _a[_i];
                    if (crit.criteria instanceof Criteria) {
                        crit.criteria.updateArrows(_this.s.criteria.length > 1);
                    }
                }
                criteria.destroy();
                _this.s.dt.draw();
                _this.s.topGroup.trigger('dtsb-redrawContents');
                return false;
            });
            criteria.dom.right
                .unbind('click')
                .on('click.dtsb', function () {
                var idx = criteria.s.index;
                var group = new Group(_this.s.dt, _this.s.opts, _this.s.topGroup, criteria.s.index, true, _this.s.depth + 1, _this.s.serverData);
                // Add the criteria that is to be moved to the new group
                group.addCriteria(criteria);
                // Update the details in the current groups criteria array
                _this.s.criteria[idx].criteria = group;
                _this.s.criteria[idx].logic = 'AND';
                _this.s.topGroup.trigger('dtsb-redrawContents');
                _this._setGroupListeners(group);
                return false;
            });
            criteria.dom.left
                .unbind('click')
                .on('click.dtsb', function () {
                _this.s.toDrop = new Criteria(_this.s.dt, _this.s.opts, _this.s.topGroup, criteria.s.index, undefined, _this.s.serverData);
                _this.s.toDrop.s = criteria.s;
                _this.s.toDrop.c = criteria.c;
                _this.s.toDrop.classes = criteria.classes;
                _this.s.toDrop.populate();
                // The dropCriteria event mutates the reference to the index so need to store it
                var index = _this.s.toDrop.s.index;
                _this.dom.container.trigger('dtsb-dropCriteria');
                criteria.s.index = index;
                _this._removeCriteria(criteria);
                // By tracking the top level group we can directly trigger a redraw on it,
                //  bubbling is also possible, but that is slow with deep levelled groups
                _this.s.topGroup.trigger('dtsb-redrawContents');
                _this.s.dt.draw();
                return false;
            });
        };
        /**
         * Set's the listeners for the group clear button
         */
        Group.prototype._setClearListener = function () {
            var _this = this;
            this.dom.clear
                .unbind('click')
                .on('click.dtsb', function () {
                if (!_this.s.isChild) {
                    _this.dom.container.trigger('dtsb-clearContents');
                    return false;
                }
                _this.destroy();
                _this.s.topGroup.trigger('dtsb-redrawContents');
                return false;
            });
        };
        /**
         * Sets listeners for sub groups of this group
         *
         * @param group The sub group that the listeners are to be set on
         */
        Group.prototype._setGroupListeners = function (group) {
            var _this = this;
            // Set listeners for the new group
            group.dom.add
                .unbind('click')
                .on('click.dtsb', function () {
                _this.setupLogic();
                _this.dom.container.trigger('dtsb-add');
                return false;
            });
            group.dom.container
                .unbind('dtsb-add')
                .on('dtsb-add.dtsb', function () {
                _this.setupLogic();
                _this.dom.container.trigger('dtsb-add');
                return false;
            });
            group.dom.container
                .unbind('dtsb-destroy')
                .on('dtsb-destroy.dtsb', function () {
                _this._removeCriteria(group, true);
                group.dom.container.remove();
                _this.setupLogic();
                return false;
            });
            group.dom.container
                .unbind('dtsb-dropCriteria')
                .on('dtsb-dropCriteria.dtsb', function () {
                var toDrop = group.s.toDrop;
                toDrop.s.index = group.s.index;
                toDrop.updateArrows(_this.s.criteria.length > 1);
                _this.addCriteria(toDrop);
                return false;
            });
            group.setListeners();
        };
        /**
         * Sets up the Group instance, setting listeners and appending elements
         */
        Group.prototype._setup = function () {
            this.setListeners();
            this.dom.add.html(this.s.dt.i18n('searchBuilder.add', this.c.i18n.add));
            this.dom.logic.children().first().html(this.c.logic === 'OR'
                ? this.s.dt.i18n('searchBuilder.logicOr', this.c.i18n.logicOr)
                : this.s.dt.i18n('searchBuilder.logicAnd', this.c.i18n.logicAnd));
            this.s.logic = this.c.logic === 'OR' ? 'OR' : 'AND';
            if (this.c.greyscale) {
                this.dom.logic.addClass(this.classes.greyscale);
            }
            this.dom.logicContainer.append(this.dom.logic).append(this.dom.clear);
            // Only append the logic button immediately if this is a sub group,
            //  otherwise it will be prepended later when adding a criteria
            if (this.s.isChild) {
                this.dom.container.append(this.dom.logicContainer);
            }
            this.dom.container.append(this.dom.add);
        };
        /**
         * Sets the listener for the logic button
         */
        Group.prototype._setLogicListener = function () {
            var _this = this;
            this.dom.logic
                .unbind('click')
                .on('click.dtsb', function () {
                _this._toggleLogic();
                _this.s.dt.draw();
                for (var _i = 0, _a = _this.s.criteria; _i < _a.length; _i++) {
                    var crit = _a[_i];
                    crit.criteria.setListeners();
                }
            });
        };
        /**
         * Toggles the logic for the group
         */
        Group.prototype._toggleLogic = function () {
            if (this.s.logic === 'OR') {
                this.s.logic = 'AND';
                this.dom.logic.children().first().html(this.s.dt.i18n('searchBuilder.logicAnd', this.c.i18n.logicAnd));
            }
            else if (this.s.logic === 'AND') {
                this.s.logic = 'OR';
                this.dom.logic.children().first().html(this.s.dt.i18n('searchBuilder.logicOr', this.c.i18n.logicOr));
            }
        };
        Group.version = '1.1.0';
        Group.classes = {
            add: 'dtsb-add',
            button: 'dtsb-button',
            clearGroup: 'dtsb-clearGroup',
            greyscale: 'dtsb-greyscale',
            group: 'dtsb-group',
            inputButton: 'dtsb-iptbtn',
            logic: 'dtsb-logic',
            logicContainer: 'dtsb-logicContainer'
        };
        Group.defaults = {
            columns: true,
            conditions: {
                'date': Criteria.dateConditions,
                'html': Criteria.stringConditions,
                'html-num': Criteria.numConditions,
                'html-num-fmt': Criteria.numFmtConditions,
                'luxon': Criteria.luxonDateConditions,
                'moment': Criteria.momentDateConditions,
                'num': Criteria.numConditions,
                'num-fmt': Criteria.numFmtConditions,
                'string': Criteria.stringConditions
            },
            depthLimit: false,
            enterSearch: false,
            filterChanged: undefined,
            greyscale: false,
            i18n: {
                add: 'Add Condition',
                button: {
                    0: 'Search Builder',
                    _: 'Search Builder (%d)'
                },
                clearAll: 'Clear All',
                condition: 'Condition',
                data: 'Data',
                "delete": '&times',
                deleteTitle: 'Delete filtering rule',
                left: '<',
                leftTitle: 'Outdent criteria',
                logicAnd: 'And',
                logicOr: 'Or',
                right: '>',
                rightTitle: 'Indent criteria',
                title: {
                    0: 'Custom Search Builder',
                    _: 'Custom Search Builder (%d)'
                },
                value: 'Value',
                valueJoiner: 'and'
            },
            logic: 'AND',
            orthogonal: {
                display: 'display',
                search: 'filter'
            },
            preDefined: false
        };
        return Group;
    }());

    var $$1;
    var dataTable$1;
    /**
     * Sets the value of jQuery for use in the file
     *
     * @param jq the instance of jQuery to be set
     */
    function setJQuery(jq) {
        $$1 = jq;
        dataTable$1 = jq.fn.DataTable;
    }
    /**
     * SearchBuilder class for DataTables.
     * Allows for complex search queries to be constructed and implemented on a DataTable
     */
    var SearchBuilder = /** @class */ (function () {
        function SearchBuilder(builderSettings, opts) {
            var _this = this;
            // Check that the required version of DataTables is included
            if (!dataTable$1 || !dataTable$1.versionCheck || !dataTable$1.versionCheck('1.10.0')) {
                throw new Error('SearchBuilder requires DataTables 1.10 or newer');
            }
            var table = new dataTable$1.Api(builderSettings);
            this.classes = $$1.extend(true, {}, SearchBuilder.classes);
            // Get options from user
            this.c = $$1.extend(true, {}, SearchBuilder.defaults, opts);
            this.dom = {
                clearAll: $$1('<button type="button">' + table.i18n('searchBuilder.clearAll', this.c.i18n.clearAll) + '</button>')
                    .addClass(this.classes.clearAll)
                    .addClass(this.classes.button)
                    .attr('type', 'button'),
                container: $$1('<div/>')
                    .addClass(this.classes.container),
                title: $$1('<div/>')
                    .addClass(this.classes.title),
                titleRow: $$1('<div/>')
                    .addClass(this.classes.titleRow),
                topGroup: undefined
            };
            this.s = {
                dt: table,
                opts: opts,
                search: undefined,
                serverData: undefined,
                topGroup: undefined
            };
            // If searchbuilder is already defined for this table then return
            if (table.settings()[0]._searchBuilder !== undefined) {
                return;
            }
            table.settings()[0]._searchBuilder = this;
            // If using SSP we want to include the previous state in the very first server call
            if (this.s.dt.page.info().serverSide) {
                this.s.dt.on('preXhr.dtsb', function (e, settings, data) {
                    var loadedState = _this.s.dt.state.loaded();
                    if (loadedState && loadedState.searchBuilder) {
                        data.searchBuilder = _this._collapseArray(loadedState.searchBuilder);
                    }
                });
                this.s.dt.on('xhr.dtsb', function (e, settings, json) {
                    if (json && json.searchBuilder && json.searchBuilder.options) {
                        _this.s.serverData = json.searchBuilder.options;
                    }
                });
            }
            // Run the remaining setup when the table is initialised
            if (this.s.dt.settings()[0]._bInitComplete) {
                this._setUp();
            }
            else {
                table.one('init.dt', function () {
                    _this._setUp();
                });
            }
            return this;
        }
        /**
         * Gets the details required to rebuild the SearchBuilder as it currently is
         */
        // eslint upset at empty object but that is what it is
        // eslint-disable-next-line @typescript-eslint/ban-types
        SearchBuilder.prototype.getDetails = function (deFormatDates) {
            if (deFormatDates === void 0) { deFormatDates = false; }
            return this.s.topGroup.getDetails(deFormatDates);
        };
        /**
         * Getter for the node of the container for the searchBuilder
         *
         * @returns JQuery<HTMLElement> the node of the container
         */
        SearchBuilder.prototype.getNode = function () {
            return this.dom.container;
        };
        /**
         * Rebuilds the SearchBuilder to a state that is provided
         *
         * @param details The details required to perform a rebuild
         */
        SearchBuilder.prototype.rebuild = function (details) {
            this.dom.clearAll.click();
            // If there are no details to rebuild then return
            if (details === undefined || details === null) {
                return this;
            }
            this.s.topGroup.s.preventRedraw = true;
            this.s.topGroup.rebuild(details);
            this.s.topGroup.s.preventRedraw = false;
            this._checkClear();
            this._updateTitle(this.s.topGroup.count());
            this.s.topGroup.redrawContents();
            this.s.dt.draw(false);
            this.s.topGroup.setListeners();
            return this;
        };
        /**
         * Applies the defaults to preDefined criteria
         *
         * @param preDef the array of criteria to be processed.
         */
        SearchBuilder.prototype._applyPreDefDefaults = function (preDef) {
            var _this = this;
            if (preDef.criteria !== undefined && preDef.logic === undefined) {
                preDef.logic = 'AND';
            }
            var _loop_1 = function (crit) {
                // Apply the defaults to any further criteria
                if (crit.criteria !== undefined) {
                    crit = this_1._applyPreDefDefaults(crit);
                }
                else {
                    this_1.s.dt.columns().every(function (index) {
                        if (_this.s.dt.settings()[0].aoColumns[index].sTitle === crit.data) {
                            crit.dataIdx = index;
                        }
                    });
                }
            };
            var this_1 = this;
            for (var _i = 0, _a = preDef.criteria; _i < _a.length; _i++) {
                var crit = _a[_i];
                _loop_1(crit);
            }
            return preDef;
        };
        /**
         * Set's up the SearchBuilder
         */
        SearchBuilder.prototype._setUp = function (loadState) {
            var _this = this;
            if (loadState === void 0) { loadState = true; }
            // Register an Api method for getting the column type
            $$1.fn.DataTable.Api.registerPlural('columns().type()', 'column().type()', function () {
                return this.iterator('column', function (settings, column) {
                    return settings.aoColumns[column].sType;
                }, 1);
            });
            // Check that DateTime is included, If not need to check if it could be used
            // eslint-disable-next-line no-extra-parens
            if (!dataTable$1.DateTime) {
                var types = this.s.dt.columns().type().toArray();
                if (types === undefined || types.includes(undefined) || types.includes(null)) {
                    types = [];
                    for (var _i = 0, _a = this.s.dt.settings()[0].aoColumns; _i < _a.length; _i++) {
                        var colInit = _a[_i];
                        types.push(colInit.searchBuilderType !== undefined ? colInit.searchBuilderType : colInit.sType);
                    }
                }
                var columnIdxs = this.s.dt.columns().toArray();
                // If the types are not yet set then draw to see if they can be retrieved then
                if (types === undefined || types.includes(undefined) || types.includes(null)) {
                    $$1.fn.dataTable.ext.oApi._fnColumnTypes(this.s.dt.settings()[0]);
                    types = this.s.dt.columns().type().toArray();
                }
                for (var i = 0; i < columnIdxs[0].length; i++) {
                    var column = columnIdxs[0][i];
                    var type = types[column];
                    if (
                    // Check if this column can be filtered
                    (this.c.columns === true ||
                        Array.isArray(this.c.columns) &&
                            this.c.columns.includes(i)) &&
                        // Check if the type is one of the restricted types
                        (type.includes('date') ||
                            type.includes('moment') ||
                            type.includes('luxon'))) {
                        alert('SearchBuilder Requires DateTime when used with dates.');
                        throw new Error('SearchBuilder requires DateTime');
                    }
                }
            }
            this.s.topGroup = new Group(this.s.dt, this.c, undefined, undefined, undefined, undefined, this.s.serverData);
            this._setClearListener();
            this.s.dt.on('stateSaveParams.dtsb', function (e, settings, data) {
                data.searchBuilder = _this.getDetails();
                if (!data.scroller) {
                    data.page = _this.s.dt.page();
                }
                else {
                    data.start = _this.s.dt.state().start;
                }
            });
            this.s.dt.on('stateLoadParams.dtsb', function (e, settings, data) {
                _this.rebuild(data.searchBuilder);
            });
            this._build();
            this.s.dt.on('preXhr.dtsb', function (e, settings, data) {
                if (_this.s.dt.page.info().serverSide) {
                    data.searchBuilder = _this._collapseArray(_this.getDetails(true));
                }
            });
            this.s.dt.on('column-reorder', function () {
                _this.rebuild(_this.getDetails());
            });
            if (loadState) {
                var loadedState = this.s.dt.state.loaded();
                // If the loaded State is not null rebuild based on it for statesave
                if (loadedState !== null && loadedState.searchBuilder !== undefined) {
                    this.s.topGroup.rebuild(loadedState.searchBuilder);
                    this.s.topGroup.dom.container.trigger('dtsb-redrawContents');
                    // If using SSP we want to restrict the amount of server calls that take place
                    //  and this information will already have been processed
                    if (!this.s.dt.page.info().serverSide) {
                        if (loadedState.page) {
                            this.s.dt.page(loadedState.page).draw('page');
                        }
                        else if (this.s.dt.scroller && loadedState.scroller) {
                            this.s.dt.scroller().scrollToRow(loadedState.scroller.topRow);
                        }
                    }
                    this.s.topGroup.setListeners();
                }
                // Otherwise load any predefined options
                else if (this.c.preDefined !== false) {
                    this.c.preDefined = this._applyPreDefDefaults(this.c.preDefined);
                    this.rebuild(this.c.preDefined);
                }
            }
            this._setEmptyListener();
            this.s.dt.state.save();
        };
        SearchBuilder.prototype._collapseArray = function (criteria) {
            if (criteria.logic === undefined) {
                if (criteria.value !== undefined) {
                    criteria.value.sort(function (a, b) {
                        if (!isNaN(+a)) {
                            a = +a;
                            b = +b;
                        }
                        if (a < b) {
                            return -1;
                        }
                        else if (b < a) {
                            return 1;
                        }
                        else {
                            return 0;
                        }
                    });
                    criteria.value1 = criteria.value[0];
                    criteria.value2 = criteria.value[1];
                }
            }
            else {
                for (var i = 0; i < criteria.criteria.length; i++) {
                    criteria.criteria[i] = this._collapseArray(criteria.criteria[i]);
                }
            }
            return criteria;
        };
        /**
         * Updates the title of the SearchBuilder
         *
         * @param count the number of filters in the SearchBuilder
         */
        SearchBuilder.prototype._updateTitle = function (count) {
            this.dom.title.html(this.s.dt.i18n('searchBuilder.title', this.c.i18n.title, count));
        };
        /**
         * Builds all of the dom elements together
         */
        SearchBuilder.prototype._build = function () {
            var _this = this;
            // Empty and setup the container
            this.dom.clearAll.remove();
            this.dom.container.empty();
            var count = this.s.topGroup.count();
            this._updateTitle(count);
            this.dom.titleRow.append(this.dom.title);
            this.dom.container.append(this.dom.titleRow);
            this.dom.topGroup = this.s.topGroup.getNode();
            this.dom.container.append(this.dom.topGroup);
            this._setRedrawListener();
            var tableNode = this.s.dt.table(0).node();
            if (!$$1.fn.dataTable.ext.search.includes(this.s.search)) {
                // Custom search function for SearchBuilder
                this.s.search = function (settings, searchData, dataIndex) {
                    if (settings.nTable !== tableNode) {
                        return true;
                    }
                    return _this.s.topGroup.search(searchData, dataIndex);
                };
                // Add SearchBuilder search function to the dataTables search array
                $$1.fn.dataTable.ext.search.push(this.s.search);
            }
            this.s.dt.on('destroy.dtsb', function () {
                _this.dom.container.remove();
                _this.dom.clearAll.remove();
                var searchIdx = $$1.fn.dataTable.ext.search.indexOf(_this.s.search);
                while (searchIdx !== -1) {
                    $$1.fn.dataTable.ext.search.splice(searchIdx, 1);
                    searchIdx = $$1.fn.dataTable.ext.search.indexOf(_this.s.search);
                }
                _this.s.dt.off('.dtsb');
                $$1(_this.s.dt.table().node()).off('.dtsb');
            });
        };
        /**
         * Checks if the clearAll button should be added or not
         */
        SearchBuilder.prototype._checkClear = function () {
            if (this.s.topGroup.s.criteria.length > 0) {
                this.dom.clearAll.insertAfter(this.dom.title);
                this._setClearListener();
            }
            else {
                this.dom.clearAll.remove();
            }
        };
        /**
         * Update the count in the title/button
         *
         * @param count Number of filters applied
         */
        SearchBuilder.prototype._filterChanged = function (count) {
            var fn = this.c.filterChanged;
            if (typeof fn === 'function') {
                fn(count, this.s.dt.i18n('searchBuilder.button', this.c.i18n.button, count));
            }
        };
        /**
         * Set the listener for the clear button
         */
        SearchBuilder.prototype._setClearListener = function () {
            var _this = this;
            this.dom.clearAll.unbind('click');
            this.dom.clearAll.on('click.dtsb', function () {
                _this.s.topGroup = new Group(_this.s.dt, _this.c, undefined, undefined, undefined, undefined, _this.s.serverData);
                _this._build();
                _this.s.dt.draw();
                _this.s.topGroup.setListeners();
                _this.dom.clearAll.remove();
                _this._setEmptyListener();
                _this._filterChanged(0);
                return false;
            });
        };
        /**
         * Set the listener for the Redraw event
         */
        SearchBuilder.prototype._setRedrawListener = function () {
            var _this = this;
            this.s.topGroup.dom.container.unbind('dtsb-redrawContents');
            this.s.topGroup.dom.container.on('dtsb-redrawContents.dtsb', function () {
                _this._checkClear();
                _this.s.topGroup.redrawContents();
                _this.s.topGroup.setupLogic();
                _this._setEmptyListener();
                var count = _this.s.topGroup.count();
                _this._updateTitle(count);
                _this._filterChanged(count);
                // If using SSP we want to restrict the amount of server calls that take place
                //  and this information will already have been processed
                if (!_this.s.dt.page.info().serverSide) {
                    _this.s.dt.draw();
                }
                _this.s.dt.state.save();
            });
            this.s.topGroup.dom.container.unbind('dtsb-redrawContents-noDraw');
            this.s.topGroup.dom.container.on('dtsb-redrawContents-noDraw.dtsb', function () {
                _this._checkClear();
                _this.s.topGroup.s.preventRedraw = true;
                _this.s.topGroup.redrawContents();
                _this.s.topGroup.s.preventRedraw = false;
                _this.s.topGroup.setupLogic();
                _this._setEmptyListener();
                var count = _this.s.topGroup.count();
                _this._updateTitle(count);
                _this._filterChanged(count);
            });
            this.s.topGroup.dom.container.unbind('dtsb-redrawLogic');
            this.s.topGroup.dom.container.on('dtsb-redrawLogic.dtsb', function () {
                _this.s.topGroup.redrawLogic();
                var count = _this.s.topGroup.count();
                _this._updateTitle(count);
                _this._filterChanged(count);
            });
            this.s.topGroup.dom.container.unbind('dtsb-add');
            this.s.topGroup.dom.container.on('dtsb-add.dtsb', function () {
                var count = _this.s.topGroup.count();
                _this._updateTitle(count);
                _this._filterChanged(count);
            });
            this.s.dt.on('postEdit.dtsb postCreate.dtsb postRemove.dtsb', function () {
                _this.s.topGroup.redrawContents();
            });
            this.s.topGroup.dom.container.unbind('dtsb-clearContents');
            this.s.topGroup.dom.container.on('dtsb-clearContents.dtsb', function () {
                _this._setUp(false);
                _this._filterChanged(0);
                _this.s.dt.draw();
            });
        };
        /**
         * Sets listeners to check whether clearAll should be added or removed
         */
        SearchBuilder.prototype._setEmptyListener = function () {
            var _this = this;
            this.s.topGroup.dom.add.on('click.dtsb', function () {
                _this._checkClear();
            });
            this.s.topGroup.dom.container.on('dtsb-destroy.dtsb', function () {
                _this.dom.clearAll.remove();
            });
        };
        SearchBuilder.version = '1.4.2';
        SearchBuilder.classes = {
            button: 'dtsb-button',
            clearAll: 'dtsb-clearAll',
            container: 'dtsb-searchBuilder',
            inputButton: 'dtsb-iptbtn',
            title: 'dtsb-title',
            titleRow: 'dtsb-titleRow'
        };
        SearchBuilder.defaults = {
            columns: true,
            conditions: {
                'date': Criteria.dateConditions,
                'html': Criteria.stringConditions,
                'html-num': Criteria.numConditions,
                'html-num-fmt': Criteria.numFmtConditions,
                'luxon': Criteria.luxonDateConditions,
                'moment': Criteria.momentDateConditions,
                'num': Criteria.numConditions,
                'num-fmt': Criteria.numFmtConditions,
                'string': Criteria.stringConditions
            },
            depthLimit: false,
            enterSearch: false,
            filterChanged: undefined,
            greyscale: false,
            i18n: {
                add: 'Add Condition',
                button: {
                    0: 'Search Builder',
                    _: 'Search Builder (%d)'
                },
                clearAll: 'Clear All',
                condition: 'Condition',
                conditions: {
                    array: {
                        contains: 'Contains',
                        empty: 'Empty',
                        equals: 'Equals',
                        not: 'Not',
                        notEmpty: 'Not Empty',
                        without: 'Without'
                    },
                    date: {
                        after: 'After',
                        before: 'Before',
                        between: 'Between',
                        empty: 'Empty',
                        equals: 'Equals',
                        not: 'Not',
                        notBetween: 'Not Between',
                        notEmpty: 'Not Empty'
                    },
                    // eslint-disable-next-line id-blacklist
                    number: {
                        between: 'Between',
                        empty: 'Empty',
                        equals: 'Equals',
                        gt: 'Greater Than',
                        gte: 'Greater Than Equal To',
                        lt: 'Less Than',
                        lte: 'Less Than Equal To',
                        not: 'Not',
                        notBetween: 'Not Between',
                        notEmpty: 'Not Empty'
                    },
                    // eslint-disable-next-line id-blacklist
                    string: {
                        contains: 'Contains',
                        empty: 'Empty',
                        endsWith: 'Ends With',
                        equals: 'Equals',
                        not: 'Not',
                        notContains: 'Does Not Contain',
                        notEmpty: 'Not Empty',
                        notEndsWith: 'Does Not End With',
                        notStartsWith: 'Does Not Start With',
                        startsWith: 'Starts With'
                    }
                },
                data: 'Data',
                "delete": '&times',
                deleteTitle: 'Delete filtering rule',
                left: '<',
                leftTitle: 'Outdent criteria',
                logicAnd: 'And',
                logicOr: 'Or',
                right: '>',
                rightTitle: 'Indent criteria',
                title: {
                    0: 'Custom Search Builder',
                    _: 'Custom Search Builder (%d)'
                },
                value: 'Value',
                valueJoiner: 'and'
            },
            logic: 'AND',
            orthogonal: {
                display: 'display',
                search: 'filter'
            },
            preDefined: false
        };
        return SearchBuilder;
    }());

    /*! SearchBuilder 1.4.2
     * ©SpryMedia Ltd - datatables.net/license/mit
     */
    setJQuery($);
    setJQuery$1($);
    setJQuery$2($);
    var dataTable = $.fn.dataTable;
    // eslint-disable-next-line no-extra-parens
    DataTable.SearchBuilder = SearchBuilder;
    // eslint-disable-next-line no-extra-parens
    dataTable.SearchBuilder = SearchBuilder;
    // eslint-disable-next-line no-extra-parens
    DataTable.Group = Group;
    // eslint-disable-next-line no-extra-parens
    dataTable.Group = Group;
    // eslint-disable-next-line no-extra-parens
    DataTable.Criteria = Criteria;
    // eslint-disable-next-line no-extra-parens
    dataTable.Criteria = Criteria;
    // eslint-disable-next-line no-extra-parens
    var apiRegister = DataTable.Api.register;
    // Set up object for plugins
    DataTable.ext.searchBuilder = {
        conditions: {}
    };
    DataTable.ext.buttons.searchBuilder = {
        action: function (e, dt, node, config) {
            this.popover(config._searchBuilder.getNode(), {
                align: 'container',
                span: 'container'
            });
            var topGroup = config._searchBuilder.s.topGroup;
            // Need to redraw the contents to calculate the correct positions for the elements
            if (topGroup !== undefined) {
                topGroup.dom.container.trigger('dtsb-redrawContents-noDraw');
            }
            if (topGroup.s.criteria.length === 0) {
                $('.' + $.fn.dataTable.Group.classes.add.replace(/ /g, '.')).click();
            }
        },
        config: {},
        init: function (dt, node, config) {
            var sb = new DataTable.SearchBuilder(dt, $.extend({
                filterChanged: function (count, text) {
                    dt.button(node).text(text);
                }
            }, config.config));
            dt.button(node).text(config.text || dt.i18n('searchBuilder.button', sb.c.i18n.button, 0));
            config._searchBuilder = sb;
        },
        text: null
    };
    apiRegister('searchBuilder.getDetails()', function (deFormatDates) {
        if (deFormatDates === void 0) { deFormatDates = false; }
        var ctx = this.context[0];
        // If SearchBuilder has not been initialised on this instance then return
        return ctx._searchBuilder ?
            ctx._searchBuilder.getDetails(deFormatDates) :
            null;
    });
    apiRegister('searchBuilder.rebuild()', function (details) {
        var ctx = this.context[0];
        // If SearchBuilder has not been initialised on this instance then return
        if (ctx._searchBuilder === undefined) {
            return null;
        }
        ctx._searchBuilder.rebuild(details);
        return this;
    });
    apiRegister('searchBuilder.container()', function () {
        var ctx = this.context[0];
        // If SearchBuilder has not been initialised on this instance then return
        return ctx._searchBuilder ?
            ctx._searchBuilder.getNode() :
            null;
    });
    /**
     * Init function for SearchBuilder
     *
     * @param settings the settings to be applied
     * @param options the options for SearchBuilder
     * @returns JQUERY<HTMLElement> Returns the node of the SearchBuilder
     */
    function _init(settings, options) {
        var api = new DataTable.Api(settings);
        var opts = options
            ? options
            : api.init().searchBuilder || DataTable.defaults.searchBuilder;
        var searchBuilder = new SearchBuilder(api, opts);
        var node = searchBuilder.getNode();
        return node;
    }
    // Attach a listener to the document which listens for DataTables initialisation
    // events so we can automatically initialise
    $(document).on('preInit.dt.dtsp', function (e, settings) {
        if (e.namespace !== 'dt') {
            return;
        }
        if (settings.oInit.searchBuilder ||
            DataTable.defaults.searchBuilder) {
            if (!settings._searchBuilder) {
                _init(settings);
            }
        }
    });
    // DataTables `dom` feature option
    DataTable.ext.feature.push({
        cFeature: 'Q',
        fnInit: _init
    });
    // DataTables 2 layout feature
    if (DataTable.ext.features) {
        DataTable.ext.features.register('searchBuilder', _init);
    }

})();


return DataTable;
}));


/*! Bootstrap 5 ui integration for DataTables' SearchBuilder
 * © SpryMedia Ltd - datatables.net/license
 */

(function( factory ){
	if ( typeof define === 'function' && define.amd ) {
		// AMD
		define( ['jquery', 'datatables.net-bs5', 'datatables.net-searchbuilder'], function ( $ ) {
			return factory( $, window, document );
		} );
	}
	else if ( typeof exports === 'object' ) {
		// CommonJS
		var jq = require('jquery');
		var cjsRequires = function (root, $) {
			if ( ! $.fn.dataTable ) {
				require('datatables.net-bs5')(root, $);
			}

			if ( ! $.fn.dataTable.SearchBuilder ) {
				require('datatables.net-searchbuilder')(root, $);
			}
		};

		if (typeof window !== 'undefined') {
			module.exports = function (root, $) {
				if ( ! root ) {
					// CommonJS environments without a window global must pass a
					// root. This will give an error otherwise
					root = window;
				}

				if ( ! $ ) {
					$ = jq( root );
				}

				cjsRequires( root, $ );
				return factory( $, root, root.document );
			};
		}
		else {
			cjsRequires( window, jq );
			module.exports = factory( jq, window, window.document );
		}
	}
	else {
		// Browser
		factory( jQuery, window, document );
	}
}(function( $, window, document, undefined ) {
'use strict';
var DataTable = $.fn.dataTable;


$.extend(true, DataTable.SearchBuilder.classes, {
    clearAll: 'btn btn-light dtsb-clearAll'
});
$.extend(true, DataTable.Group.classes, {
    add: 'btn btn-light dtsb-add',
    clearGroup: 'btn btn-light dtsb-clearGroup',
    logic: 'btn btn-light dtsb-logic'
});
$.extend(true, DataTable.Criteria.classes, {
    condition: 'form-select dtsb-condition',
    data: 'dtsb-data form-select',
    "delete": 'btn btn-light dtsb-delete',
    input: 'form-control dtsb-input',
    left: 'btn btn-light dtsb-left',
    right: 'btn btn-light dtsb-right',
    select: 'form-select',
    value: 'dtsb-value'
});


return DataTable;
}));


/*! SearchPanes 2.1.2
 * © SpryMedia Ltd - datatables.net/license
 */

(function( factory ){
	if ( typeof define === 'function' && define.amd ) {
		// AMD
		define( ['jquery', 'datatables.net'], function ( $ ) {
			return factory( $, window, document );
		} );
	}
	else if ( typeof exports === 'object' ) {
		// CommonJS
		var jq = require('jquery');
		var cjsRequires = function (root, $) {
			if ( ! $.fn.dataTable ) {
				require('datatables.net')(root, $);
			}
		};

		if (typeof window !== 'undefined') {
			module.exports = function (root, $) {
				if ( ! root ) {
					// CommonJS environments without a window global must pass a
					// root. This will give an error otherwise
					root = window;
				}

				if ( ! $ ) {
					$ = jq( root );
				}

				cjsRequires( root, $ );
				return factory( $, root, root.document );
			};
		}
		else {
			cjsRequires( window, jq );
			module.exports = factory( jq, window, window.document );
		}
	}
	else {
		// Browser
		factory( jQuery, window, document );
	}
}(function( $, window, document, undefined ) {
'use strict';
var DataTable = $.fn.dataTable;


(function () {
    'use strict';

    var $$5;
    var dataTable$2;
    function setJQuery$4(jq) {
        $$5 = jq;
        dataTable$2 = jq.fn.dataTable;
    }
    var SearchPane = /** @class */ (function () {
        /**
         * Creates the panes, sets up the search function
         *
         * @param paneSettings The settings for the searchPanes
         * @param opts The options for the default features
         * @param index the index of the column for this pane
         * @param panesContainer The overall container for SearchPanes that this pane will be attached to
         * @param panes The custom pane settings if this is a custom pane
         * @returns {object} the pane that has been created, including the table and the index of the pane
         */
        function SearchPane(paneSettings, opts, index, panesContainer, panes) {
            var _this = this;
            if (panes === void 0) { panes = null; }
            // Check that the required version of DataTables is included
            if (!dataTable$2 || !dataTable$2.versionCheck || !dataTable$2.versionCheck('1.10.0')) {
                throw new Error('SearchPane requires DataTables 1.10 or newer');
            }
            // Check that Select is included
            // eslint-disable-next-line no-extra-parens
            if (!dataTable$2.select) {
                throw new Error('SearchPane requires Select');
            }
            var table = new dataTable$2.Api(paneSettings);
            this.classes = $$5.extend(true, {}, SearchPane.classes);
            // Get options from user
            this.c = $$5.extend(true, {}, SearchPane.defaults, opts, panes);
            if (opts && opts.hideCount && opts.viewCount === undefined) {
                this.c.viewCount = !this.c.hideCount;
            }
            var rowLength = table.columns().eq(0).toArray().length;
            this.s = {
                colExists: index < rowLength,
                colOpts: undefined,
                customPaneSettings: panes,
                displayed: false,
                dt: table,
                dtPane: undefined,
                firstSet: true,
                index: index,
                indexes: [],
                listSet: false,
                name: undefined,
                rowData: {
                    arrayFilter: [],
                    arrayOriginal: [],
                    bins: {},
                    binsOriginal: {},
                    filterMap: new Map(),
                    totalOptions: 0
                },
                scrollTop: 0,
                searchFunction: undefined,
                selections: [],
                serverSelect: [],
                serverSelecting: false,
                tableLength: null,
                updating: false
            };
            this.s.colOpts = this.s.colExists ? this._getOptions() : this._getBonusOptions();
            this.dom = {
                buttonGroup: $$5('<div/>').addClass(this.classes.buttonGroup),
                clear: $$5('<button type="button">&#215;</button>')
                    .attr('disabled', 'true')
                    .addClass(this.classes.disabledButton)
                    .addClass(this.classes.paneButton)
                    .addClass(this.classes.clearButton)
                    .html(this.s.dt.i18n('searchPanes.clearPane', this.c.i18n.clearPane)),
                collapseButton: $$5('<button type="button"><span class="' + this.classes.caret + '">&#x5e;</span></button>')
                    .addClass(this.classes.paneButton)
                    .addClass(this.classes.collapseButton),
                container: $$5('<div/>')
                    .addClass(this.classes.container)
                    .addClass(this.s.colOpts.className)
                    .addClass(this.classes.layout +
                    (parseInt(this.c.layout.split('-')[1], 10) < 10 ?
                        this.c.layout :
                        this.c.layout.split('-')[0] + '-9'))
                    .addClass(this.s.customPaneSettings && this.s.customPaneSettings.className
                    ? this.s.customPaneSettings.className
                    : ''),
                countButton: $$5('<button type="button"></button>')
                    .addClass(this.classes.paneButton)
                    .addClass(this.classes.countButton),
                dtP: $$5('<table><thead><tr><th>' +
                    (this.s.colExists
                        ? $$5(this.s.dt.column(this.s.index).header()).text()
                        : this.s.customPaneSettings.header || 'Custom Pane') + '</th><th/></tr></thead></table>'),
                lower: $$5('<div/>').addClass(this.classes.subRow2).addClass(this.classes.narrowButton),
                nameButton: $$5('<button type="button"></button>')
                    .addClass(this.classes.paneButton)
                    .addClass(this.classes.nameButton),
                panesContainer: panesContainer,
                searchBox: $$5('<input/>').addClass(this.classes.paneInputButton).addClass(this.classes.search),
                searchButton: $$5('<button type = "button"/>')
                    .addClass(this.classes.searchIcon)
                    .addClass(this.classes.paneButton),
                searchCont: $$5('<div/>').addClass(this.classes.searchCont),
                searchLabelCont: $$5('<div/>').addClass(this.classes.searchLabelCont),
                topRow: $$5('<div/>').addClass(this.classes.topRow),
                upper: $$5('<div/>').addClass(this.classes.subRow1).addClass(this.classes.narrowSearch)
            };
            // Set the value of name incase ordering is desired
            if (this.s.colOpts.name) {
                this.s.name = this.s.colOpts.name;
            }
            else if (this.s.customPaneSettings && this.s.customPaneSettings.name) {
                this.s.name = this.s.customPaneSettings.name;
            }
            else {
                this.s.name = this.s.colExists ?
                    $$5(this.s.dt.column(this.s.index).header()).text() :
                    this.s.customPaneSettings.header || 'Custom Pane';
            }
            var tableNode = this.s.dt.table(0).node();
            // Custom search function for table
            this.s.searchFunction = function (settings, searchData, dataIndex) {
                // If no data has been selected then show all
                if (_this.s.selections.length === 0) {
                    return true;
                }
                if (settings.nTable !== tableNode) {
                    return true;
                }
                var filter = null;
                if (_this.s.colExists) {
                    // Get the current filtered data
                    filter = searchData[_this.s.index];
                    if (_this.s.colOpts.orthogonal.filter !== 'filter') {
                        // get the filter value from the map
                        filter = _this.s.rowData.filterMap.get(dataIndex);
                        if (filter instanceof $$5.fn.dataTable.Api) {
                            // eslint-disable-next-line no-extra-parens
                            filter = filter.toArray();
                        }
                    }
                }
                return _this._search(filter, dataIndex);
            };
            $$5.fn.dataTable.ext.search.push(this.s.searchFunction);
            // If the clear button for this pane is clicked clear the selections
            if (this.c.clear) {
                this.dom.clear.on('click.dtsp', function () {
                    var searches = _this.dom.container.find('.' + _this.classes.search.replace(/\s+/g, '.'));
                    searches.each(function () {
                        $$5(this).val('').trigger('input');
                    });
                    _this.clearPane();
                });
            }
            // Sometimes the top row of the panes containing the search box and ordering buttons appears
            //  weird if the width of the panes is lower than expected, this fixes the design.
            // Equally this may occur when the table is resized.
            this.s.dt.on('draw.dtsp', function () { return _this.adjustTopRow(); });
            this.s.dt.on('buttons-action.dtsp', function () { return _this.adjustTopRow(); });
            // When column-reorder is present and the columns are moved, it is necessary to
            //  reassign all of the panes indexes to the new index of the column.
            this.s.dt.on('column-reorder.dtsp', function (e, settings, details) {
                _this.s.index = details.mapping[_this.s.index];
            });
            return this;
        }
        /**
         * Adds a row to the panes table
         *
         * @param display the value to be displayed to the user
         * @param filter the value to be filtered on when searchpanes is implemented
         * @param shown the number of rows in the table that are currently visible matching this criteria
         * @param total the total number of rows in the table that match this criteria
         * @param sort the value to be sorted in the pane table
         * @param type the value of which the type is to be derived from
         */
        SearchPane.prototype.addRow = function (display, filter, sort, type, className, total, shown) {
            if (!total) {
                total = this.s.rowData.bins[filter] ?
                    this.s.rowData.bins[filter] :
                    0;
            }
            if (!shown) {
                shown = this._getShown(filter);
            }
            var index;
            for (var _i = 0, _a = this.s.indexes; _i < _a.length; _i++) {
                var entry = _a[_i];
                if (entry.filter === filter) {
                    index = entry.index;
                }
            }
            if (index === undefined) {
                index = this.s.indexes.length;
                this.s.indexes.push({ filter: filter, index: index });
            }
            return this.s.dtPane.row.add({
                className: className,
                display: display !== '' ?
                    display :
                    this.emptyMessage(),
                filter: filter,
                index: index,
                shown: shown,
                sort: sort,
                total: total,
                type: type
            });
        };
        /**
         * Adjusts the layout of the top row when the screen is resized
         */
        SearchPane.prototype.adjustTopRow = function () {
            var subContainers = this.dom.container.find('.' + this.classes.subRowsContainer.replace(/\s+/g, '.'));
            var subRow1 = this.dom.container.find('.' + this.classes.subRow1.replace(/\s+/g, '.'));
            var subRow2 = this.dom.container.find('.' + this.classes.subRow2.replace(/\s+/g, '.'));
            var topRow = this.dom.container.find('.' + this.classes.topRow.replace(/\s+/g, '.'));
            // If the width is 0 then it is safe to assume that the pane has not yet been displayed.
            //  Even if it has, if the width is 0 it won't make a difference if it has the narrow class or not
            if (($$5(subContainers[0]).width() < 252 || $$5(topRow[0]).width() < 252) && $$5(subContainers[0]).width() !== 0) {
                $$5(subContainers[0]).addClass(this.classes.narrow);
                $$5(subRow1[0]).addClass(this.classes.narrowSub).removeClass(this.classes.narrowSearch);
                $$5(subRow2[0]).addClass(this.classes.narrowSub).removeClass(this.classes.narrowButton);
            }
            else {
                $$5(subContainers[0]).removeClass(this.classes.narrow);
                $$5(subRow1[0]).removeClass(this.classes.narrowSub).addClass(this.classes.narrowSearch);
                $$5(subRow2[0]).removeClass(this.classes.narrowSub).addClass(this.classes.narrowButton);
            }
        };
        /**
         * In the case of a rebuild there is potential for new data to have been included or removed
         * so all of the rowData must be reset as a precaution.
         */
        SearchPane.prototype.clearData = function () {
            this.s.rowData = {
                arrayFilter: [],
                arrayOriginal: [],
                bins: {},
                binsOriginal: {},
                filterMap: new Map(),
                totalOptions: 0
            };
        };
        /**
         * Clear the selections in the pane
         */
        SearchPane.prototype.clearPane = function () {
            // Deselect all rows which are selected and update the table and filter count.
            this.s.dtPane.rows({ selected: true }).deselect();
            this.updateTable();
            return this;
        };
        /**
         * Collapses the pane so that only the header is displayed
         */
        SearchPane.prototype.collapse = function () {
            var _this = this;
            if (!this.s.displayed ||
                (
                // If collapsing is disabled globally, and not enabled specifically for this column
                !this.c.collapse && this.s.colOpts.collapse !== true ||
                    // OR, collapsing could be enabled globally and this column specifically
                    // is not to be collapsed.
                    // We can't just take !this.s.colOpts.collapse here as if it is undefined
                    // then the global should be used
                    this.s.colOpts.collapse === false)) {
                return;
            }
            $$5(this.s.dtPane.table().container()).addClass(this.classes.hidden);
            this.dom.topRow.addClass(this.classes.bordered);
            this.dom.nameButton.addClass(this.classes.disabledButton);
            this.dom.countButton.addClass(this.classes.disabledButton);
            this.dom.searchButton.addClass(this.classes.disabledButton);
            this.dom.collapseButton.addClass(this.classes.rotated);
            this.dom.topRow.one('click.dtsp', function () { return _this.show(); });
            this.dom.topRow.trigger('collapse.dtsps');
        };
        /**
         * Strips all of the SearchPanes elements from the document and turns all of the listeners for the buttons off
         */
        SearchPane.prototype.destroy = function () {
            if (this.s.dtPane) {
                this.s.dtPane.off('.dtsp');
            }
            this.s.dt.off('.dtsp');
            this.dom.clear.off('.dtsp');
            this.dom.nameButton.off('.dtsp');
            this.dom.countButton.off('.dtsp');
            this.dom.searchButton.off('.dtsp');
            this.dom.collapseButton.off('.dtsp');
            $$5(this.s.dt.table().node()).off('.dtsp');
            this.dom.container.detach();
            var searchIdx = $$5.fn.dataTable.ext.search.indexOf(this.s.searchFunction);
            while (searchIdx !== -1) {
                $$5.fn.dataTable.ext.search.splice(searchIdx, 1);
                searchIdx = $$5.fn.dataTable.ext.search.indexOf(this.s.searchFunction);
            }
            // If the datatables have been defined for the panes then also destroy these
            if (this.s.dtPane) {
                this.s.dtPane.destroy();
            }
            this.s.listSet = false;
        };
        /**
         * Getting the legacy message is a little complex due a legacy parameter
         */
        SearchPane.prototype.emptyMessage = function () {
            var def = this.c.i18n.emptyMessage;
            // Legacy parameter support
            if (this.c.emptyMessage) {
                def = this.c.emptyMessage;
            }
            // Override per column
            if (this.s.colOpts.emptyMessage !== false && this.s.colOpts.emptyMessage !== null) {
                def = this.s.colOpts.emptyMessage;
            }
            return this.s.dt.i18n('searchPanes.emptyMessage', def);
        };
        /**
         * Updates the number of filters that have been applied in the title
         */
        SearchPane.prototype.getPaneCount = function () {
            return this.s.dtPane ?
                this.s.dtPane.rows({ selected: true }).data().toArray().length :
                0;
        };
        /**
         * Rebuilds the panes from the start having deleted the old ones
         *
         * @param? dataIn data to be used in buildPane
         * @param? maintainSelection Whether the current selections are to be maintained over rebuild
         */
        SearchPane.prototype.rebuildPane = function (dataIn, maintainSelection) {
            if (dataIn === void 0) { dataIn = null; }
            if (maintainSelection === void 0) { maintainSelection = false; }
            this.clearData();
            var selectedRows = [];
            this.s.serverSelect = [];
            var prevEl = null;
            // When rebuilding strip all of the HTML Elements out of the container and start from scratch
            if (this.s.dtPane) {
                if (maintainSelection) {
                    if (!this.s.dt.page.info().serverSide) {
                        selectedRows = this.s.dtPane.rows({ selected: true }).data().toArray();
                    }
                    else {
                        this.s.serverSelect = this.s.dtPane.rows({ selected: true }).data().toArray();
                    }
                }
                this.s.dtPane.clear().destroy();
                prevEl = this.dom.container.prev();
                this.destroy();
                this.s.dtPane = undefined;
                $$5.fn.dataTable.ext.search.push(this.s.searchFunction);
            }
            this.dom.container.removeClass(this.classes.hidden);
            this.s.displayed = false;
            this._buildPane(!this.s.dt.page.info().serverSide ?
                selectedRows :
                this.s.serverSelect, dataIn, prevEl);
            return this;
        };
        /**
         * Resizes the pane based on the layout that is passed in
         *
         * @param layout the layout to be applied to this pane
         */
        SearchPane.prototype.resize = function (layout) {
            this.c.layout = layout;
            this.dom.container
                .removeClass()
                .addClass(this.classes.show)
                .addClass(this.classes.container)
                .addClass(this.s.colOpts.className)
                .addClass(this.classes.layout +
                (parseInt(layout.split('-')[1], 10) < 10 ?
                    layout :
                    layout.split('-')[0] + '-9'))
                .addClass(this.s.customPaneSettings !== null && this.s.customPaneSettings.className
                ? this.s.customPaneSettings.className
                : '');
            this.adjustTopRow();
        };
        /**
         * Sets the listeners for the pane.
         *
         * Having it in it's own function makes it easier to only set them once
         */
        SearchPane.prototype.setListeners = function () {
            var _this = this;
            if (!this.s.dtPane) {
                return;
            }
            // When an item is selected on the pane, add these to the array which holds selected items.
            // Custom search will perform.
            this.s.dtPane.off('select.dtsp').on('select.dtsp', function () {
                clearTimeout(_this.s.deselectTimeout);
                _this._updateSelection(!_this.s.updating);
                _this.dom.clear.removeClass(_this.classes.disabledButton).removeAttr('disabled');
            });
            // When an item is deselected on the pane, re add the currently selected items to the array
            // which holds selected items. Custom search will be performed.
            this.s.dtPane.off('deselect.dtsp').on('deselect.dtsp', function () {
                _this.s.deselectTimeout = setTimeout(function () {
                    _this._updateSelection(true);
                    if (_this.s.dtPane.rows({ selected: true }).data().toArray().length === 0) {
                        _this.dom.clear.addClass(_this.classes.disabledButton).attr('disabled', 'true');
                    }
                }, 50);
            });
            // If we attempty to turn off this event then it will ruin behaviour in other panes
            //  so need to make sure that it is only done once
            if (this.s.firstSet) {
                this.s.firstSet = false;
                // When saving the state store all of the selected rows for preselection next time around
                this.s.dt.on('stateSaveParams.dtsp', function (e, settings, data) {
                    // If the data being passed in is empty then state clear must have occured
                    // so clear the panes state as well
                    if ($$5.isEmptyObject(data)) {
                        _this.s.dtPane.state.clear();
                        return;
                    }
                    var bins;
                    var order;
                    var selected = [];
                    var collapsed;
                    var searchTerm;
                    var arrayFilter;
                    // Get all of the data needed for the state save from the pane
                    if (_this.s.dtPane) {
                        selected = _this.s.dtPane
                            .rows({ selected: true })
                            .data()
                            .map(function (item) { return item.filter.toString(); })
                            .toArray();
                        searchTerm = _this.dom.searchBox.val();
                        order = _this.s.dtPane.order();
                        bins = _this.s.rowData.binsOriginal;
                        arrayFilter = _this.s.rowData.arrayOriginal;
                        collapsed = _this.dom.collapseButton.hasClass(_this.classes.rotated);
                    }
                    if (data.searchPanes === undefined) {
                        data.searchPanes = {};
                    }
                    if (data.searchPanes.panes === undefined) {
                        data.searchPanes.panes = [];
                    }
                    for (var i = 0; i < data.searchPanes.panes.length; i++) {
                        if (data.searchPanes.panes[i].id === _this.s.index) {
                            data.searchPanes.panes.splice(i, 1);
                            i--;
                        }
                    }
                    // Add the panes data to the state object
                    data.searchPanes.panes.push({
                        arrayFilter: arrayFilter,
                        bins: bins,
                        collapsed: collapsed,
                        id: _this.s.index,
                        order: order,
                        searchTerm: searchTerm,
                        selected: selected
                    });
                });
            }
            this.s.dtPane.off('user-select.dtsp').on('user-select.dtsp', function (e, _dt, type, cell, originalEvent) {
                originalEvent.stopPropagation();
            });
            this.s.dtPane.off('draw.dtsp').on('draw.dtsp', function () { return _this.adjustTopRow(); });
            // When the button to order by the name of the options is clicked then
            //  change the ordering to whatever it isn't currently
            this.dom.nameButton.off('click.dtsp').on('click.dtsp', function () {
                var currentOrder = _this.s.dtPane.order()[0][1];
                _this.s.dtPane.order([0, currentOrder === 'asc' ? 'desc' : 'asc']).draw();
                // This state save is required so that the ordering of the panes is maintained
                _this.s.dt.state.save();
            });
            // When the button to order by the number of entries in the column is clicked then
            //  change the ordering to whatever it isn't currently
            this.dom.countButton.off('click.dtsp').on('click.dtsp', function () {
                var currentOrder = _this.s.dtPane.order()[0][1];
                _this.s.dtPane.order([1, currentOrder === 'asc' ? 'desc' : 'asc']).draw();
                // This state save is required so that the ordering of the panes is maintained
                _this.s.dt.state.save();
            });
            // When the button to order by the number of entries in the column is clicked then
            //  change the ordering to whatever it isn't currently
            this.dom.collapseButton.off('click.dtsp').on('click.dtsp', function (e) {
                e.stopPropagation();
                var container = $$5(_this.s.dtPane.table().container());
                // Toggle the classes
                container.toggleClass(_this.classes.hidden);
                _this.dom.topRow.toggleClass(_this.classes.bordered);
                _this.dom.nameButton.toggleClass(_this.classes.disabledButton);
                _this.dom.countButton.toggleClass(_this.classes.disabledButton);
                _this.dom.searchButton.toggleClass(_this.classes.disabledButton);
                _this.dom.collapseButton.toggleClass(_this.classes.rotated);
                if (container.hasClass(_this.classes.hidden)) {
                    _this.dom.topRow.on('click.dtsp', function () { return _this.dom.collapseButton.click(); });
                }
                else {
                    _this.dom.topRow.off('click.dtsp');
                }
                _this.s.dt.state.save();
                _this.dom.topRow.trigger('collapse.dtsps');
            });
            // When the clear button is clicked reset the pane
            this.dom.clear.off('click.dtsp').on('click.dtsp', function () {
                var searches = _this.dom.container.find('.' + _this.classes.search.replace(/ /g, '.'));
                searches.each(function () {
                    // set the value of the search box to be an empty string and then search on that, effectively reseting
                    $$5(this).val('').trigger('input');
                });
                _this.clearPane();
            });
            // When the search button is clicked then draw focus to the search box
            this.dom.searchButton.off('click.dtsp').on('click.dtsp', function () { return _this.dom.searchBox.focus(); });
            // When a character is inputted into the searchbox search the pane for matching values.
            // Doing it this way means that no button has to be clicked to trigger a search, it is done asynchronously
            this.dom.searchBox.off('click.dtsp').on('input.dtsp', function () {
                var searchval = _this.dom.searchBox.val();
                _this.s.dtPane.search(searchval).draw();
                if (typeof searchval === 'string' &&
                    (searchval.length > 0 ||
                        searchval.length === 0 && _this.s.dtPane.rows({ selected: true }).data().toArray().length > 0)) {
                    _this.dom.clear.removeClass(_this.classes.disabledButton).removeAttr('disabled');
                }
                else {
                    _this.dom.clear.addClass(_this.classes.disabledButton).attr('disabled', 'true');
                }
                // This state save is required so that the searching on the panes is maintained
                _this.s.dt.state.save();
            });
            this.s.dtPane.select.style(this.s.colOpts.dtOpts && this.s.colOpts.dtOpts.select && this.s.colOpts.dtOpts.select.style
                ? this.s.colOpts.dtOpts.select.style
                : this.c.dtOpts && this.c.dtOpts.select && this.c.dtOpts.select.style
                    ? this.c.dtOpts.select.style
                    : 'os');
        };
        /**
         * Populates the SearchPane based off of the data that has been recieved from the server
         *
         * This method is overriden by SearchPaneST
         *
         * @param dataIn The data that has been sent from the server
         */
        SearchPane.prototype._serverPopulate = function (dataIn) {
            if (dataIn.tableLength) {
                this.s.tableLength = dataIn.tableLength;
                this.s.rowData.totalOptions = this.s.tableLength;
            }
            else if (this.s.tableLength === null || this.s.dt.rows()[0].length > this.s.tableLength) {
                this.s.tableLength = this.s.dt.rows()[0].length;
                this.s.rowData.totalOptions = this.s.tableLength;
            }
            var colTitle = this.s.dt.column(this.s.index).dataSrc();
            // If there is SP data for this column add it to the data array and bin
            if (dataIn.searchPanes.options[colTitle]) {
                for (var _i = 0, _a = dataIn.searchPanes.options[colTitle]; _i < _a.length; _i++) {
                    var dataPoint = _a[_i];
                    this.s.rowData.arrayFilter.push({
                        display: dataPoint.label,
                        filter: dataPoint.value,
                        sort: dataPoint.label,
                        type: dataPoint.label
                    });
                    this.s.rowData.bins[dataPoint.value] = dataPoint.total;
                }
            }
            var binLength = Object.keys(this.s.rowData.bins).length;
            var uniqueRatio = this._uniqueRatio(binLength, this.s.tableLength);
            // Don't show the pane if there isnt enough variance in the data, or there is only 1 entry for that pane
            if (this.s.displayed === false &&
                ((this.s.colOpts.show === undefined && this.s.colOpts.threshold === null ?
                    uniqueRatio > this.c.threshold :
                    uniqueRatio > this.s.colOpts.threshold) ||
                    this.s.colOpts.show !== true && binLength <= 1)) {
                this.dom.container.addClass(this.classes.hidden);
                this.s.displayed = false;
                return;
            }
            // Store the original data
            this.s.rowData.arrayOriginal = this.s.rowData.arrayFilter;
            this.s.rowData.binsOriginal = this.s.rowData.bins;
            // Flag this pane as being displayed
            this.s.displayed = true;
        };
        /**
         * Expands the pane from the collapsed state
         */
        SearchPane.prototype.show = function () {
            if (!this.s.displayed) {
                return;
            }
            this.dom.topRow.removeClass(this.classes.bordered);
            this.dom.nameButton.removeClass(this.classes.disabledButton);
            this.dom.countButton.removeClass(this.classes.disabledButton);
            this.dom.searchButton.removeClass(this.classes.disabledButton);
            this.dom.collapseButton.removeClass(this.classes.rotated);
            $$5(this.s.dtPane.table().container()).removeClass(this.classes.hidden);
            this.dom.topRow.trigger('collapse.dtsps');
        };
        /**
         * Finds the ratio of the number of different options in the table to the number of rows
         *
         * @param bins the number of different options in the table
         * @param rowCount the total number of rows in the table
         * @returns {number} returns the ratio
         */
        SearchPane.prototype._uniqueRatio = function (bins, rowCount) {
            if (rowCount > 0 &&
                (this.s.rowData.totalOptions > 0 && !this.s.dt.page.info().serverSide ||
                    this.s.dt.page.info().serverSide && this.s.tableLength > 0)) {
                return bins / this.s.rowData.totalOptions;
            }
            return 1;
        };
        /**
         * Updates the panes if one of the options to do so has been set to true
         * rather than the filtered message when using viewTotal.
         */
        SearchPane.prototype.updateTable = function () {
            var selectedRows = this.s.dtPane.rows({ selected: true }).data().toArray().map(function (el) { return el.filter; });
            this.s.selections = selectedRows;
            this._searchExtras();
        };
        /**
         * Adds the custom options to the pane
         *
         * @returns {Array} Returns the array of rows which have been added to the pane
         */
        SearchPane.prototype._getComparisonRows = function () {
            // Find the appropriate options depending on whether this is a pane for a specific column or a custom pane
            var options = this.s.colOpts.options
                ? this.s.colOpts.options
                : this.s.customPaneSettings && this.s.customPaneSettings.options
                    ? this.s.customPaneSettings.options
                    : undefined;
            if (options === undefined) {
                return;
            }
            var allRows = this.s.dt.rows();
            var tableValsTotal = allRows.data().toArray();
            var rows = [];
            // Clear all of the other rows from the pane, only custom options are to be displayed when they are defined
            this.s.dtPane.clear();
            this.s.indexes = [];
            for (var _i = 0, options_1 = options; _i < options_1.length; _i++) {
                var comp = options_1[_i];
                // Initialise the object which is to be placed in the row
                var insert = comp.label !== '' ?
                    comp.label :
                    this.emptyMessage();
                var comparisonObj = {
                    className: comp.className,
                    display: insert,
                    filter: typeof comp.value === 'function' ? comp.value : [],
                    sort: insert,
                    total: 0,
                    type: insert
                };
                // If a custom function is in place
                if (typeof comp.value === 'function') {
                    // Count the number of times the function evaluates to true for the original data in the Table
                    for (var i = 0; i < tableValsTotal.length; i++) {
                        if (comp.value.call(this.s.dt, tableValsTotal[i], allRows[0][i])) {
                            comparisonObj.total++;
                        }
                    }
                    // Update the comparisonObj
                    if (typeof comparisonObj.filter !== 'function') {
                        comparisonObj.filter.push(comp.filter);
                    }
                }
                rows.push(this.addRow(comparisonObj.display, comparisonObj.filter, comparisonObj.sort, comparisonObj.type, comparisonObj.className, comparisonObj.total));
            }
            return rows;
        };
        SearchPane.prototype._getMessage = function (row) {
            return this.s.dt.i18n('searchPanes.count', this.c.i18n.count).replace(/{total}/g, row.total);
        };
        /**
         * Overridden in SearchPaneViewTotal and SearchPaneCascade to get the number of times a specific value is shown
         *
         * Here it is blanked so that it takes no action
         *
         * @param filter The filter value
         * @returns undefined
         */
        // eslint-disable-next-line @typescript-eslint/no-unused-vars
        SearchPane.prototype._getShown = function (filter) {
            return undefined;
        };
        /**
         * Get's the pane config appropriate to this class
         *
         * @returns The config needed to create a pane of this type
         */
        SearchPane.prototype._getPaneConfig = function () {
            var _this = this;
            // eslint-disable-next-line no-extra-parens
            var haveScroller = dataTable$2.Scroller;
            var langOpts = this.s.dt.settings()[0].oLanguage;
            langOpts.url = undefined;
            langOpts.sUrl = undefined;
            return {
                columnDefs: [
                    {
                        className: 'dtsp-nameColumn',
                        data: 'display',
                        render: function (data, type, row) {
                            if (type === 'sort') {
                                return row.sort;
                            }
                            else if (type === 'type') {
                                return row.type;
                            }
                            var message = _this._getMessage(row);
                            // We are displaying the count in the same columne as the name of the search option.
                            // This is so that there is not need to call columns.adjust()
                            //  which in turn speeds up the code
                            var pill = '<span class="' + _this.classes.pill + '">' + message + '</span>';
                            if (!_this.c.viewCount || !_this.s.colOpts.viewCount) {
                                pill = '';
                            }
                            if (type === 'filter') {
                                return typeof data === 'string' && data.match(/<[^>]*>/) !== null ?
                                    data.replace(/<[^>]*>/g, '') :
                                    data;
                            }
                            return '<div class="' + _this.classes.nameCont + '"><span title="' +
                                (typeof data === 'string' && data.match(/<[^>]*>/) !== null ?
                                    data.replace(/<[^>]*>/g, '') :
                                    data) +
                                '" class="' + _this.classes.name + '">' +
                                data + '</span>' +
                                pill + '</div>';
                        },
                        targets: 0,
                        // Accessing the private datatables property to set type based on the original table.
                        // This is null if not defined by the user, meaning that automatic type detection
                        //  would take place
                        type: this.s.dt.settings()[0].aoColumns[this.s.index] ?
                            this.s.dt.settings()[0].aoColumns[this.s.index]._sManualType :
                            null
                    },
                    {
                        className: 'dtsp-countColumn ' + this.classes.badgePill,
                        data: 'total',
                        searchable: false,
                        targets: 1,
                        visible: false
                    }
                ],
                deferRender: true,
                dom: 't',
                info: false,
                language: langOpts,
                paging: haveScroller ? true : false,
                scrollX: false,
                scrollY: '200px',
                scroller: haveScroller ? true : false,
                select: true,
                stateSave: this.s.dt.settings()[0].oFeatures.bStateSave ? true : false
            };
        };
        /**
         * This method allows for changes to the panes and table to be made when a selection or a deselection occurs
         */
        SearchPane.prototype._makeSelection = function () {
            this.updateTable();
            this.s.updating = true;
            this.s.dt.draw(false);
            this.s.updating = false;
        };
        /**
         * Populates an array with all of the data for the table
         *
         * @param rowIdx The current row index to be compared
         * @param arrayFilter The array that is to be populated with row Details
         * @param settings The DataTable settings object
         * @param bins The bins object that is to be populated with the row counts
         */
        SearchPane.prototype._populatePaneArray = function (rowIdx, arrayFilter, settings, bins) {
            if (bins === void 0) { bins = this.s.rowData.bins; }
            // Retrieve the rendered data from the cell using the fnGetCellData function
            // rather than the cell().render API method for optimisation
            if (typeof this.s.colOpts.orthogonal === 'string') {
                var rendered = settings.oApi._fnGetCellData(settings, rowIdx, this.s.index, this.s.colOpts.orthogonal);
                this.s.rowData.filterMap.set(rowIdx, rendered);
                this._addOption(rendered, rendered, rendered, rendered, arrayFilter, bins);
                this.s.rowData.totalOptions++;
            }
            else {
                var filter = settings.oApi._fnGetCellData(settings, rowIdx, this.s.index, this.s.colOpts.orthogonal.search);
                // Null and empty string are to be considered the same value
                if (filter === null) {
                    filter = '';
                }
                if (typeof filter === 'string') {
                    filter = filter.replace(/<[^>]*>/g, '');
                }
                this.s.rowData.filterMap.set(rowIdx, filter);
                if (!bins[filter]) {
                    bins[filter] = 1;
                    this._addOption(filter, settings.oApi._fnGetCellData(settings, rowIdx, this.s.index, this.s.colOpts.orthogonal.display), settings.oApi._fnGetCellData(settings, rowIdx, this.s.index, this.s.colOpts.orthogonal.sort), settings.oApi._fnGetCellData(settings, rowIdx, this.s.index, this.s.colOpts.orthogonal.type), arrayFilter, bins);
                    this.s.rowData.totalOptions++;
                }
                else {
                    bins[filter]++;
                    this.s.rowData.totalOptions++;
                }
            }
        };
        /**
         * Reloads all of the previous selects into the panes
         *
         * @param loadedFilter The loaded filters from a previous state
         */
        SearchPane.prototype._reloadSelect = function (loadedFilter) {
            // If the state was not saved don't selected any
            if (loadedFilter === undefined) {
                return;
            }
            var idx;
            // For each pane, check that the loadedFilter list exists and is not null,
            // find the id of each search item and set it to be selected.
            for (var i = 0; i < loadedFilter.searchPanes.panes.length; i++) {
                if (loadedFilter.searchPanes.panes[i].id === this.s.index) {
                    idx = i;
                    break;
                }
            }
            if (idx) {
                var table = this.s.dtPane;
                var rows = table.rows({ order: 'index' }).data().map(function (item) { return item.filter !== null ?
                    item.filter.toString() :
                    null; }).toArray();
                for (var _i = 0, _a = loadedFilter.searchPanes.panes[idx].selected; _i < _a.length; _i++) {
                    var filter = _a[_i];
                    var id = -1;
                    if (filter !== null) {
                        id = rows.indexOf(filter.toString());
                    }
                    if (id > -1) {
                        this.s.serverSelecting = true;
                        table.row(id).select();
                        this.s.serverSelecting = false;
                    }
                }
            }
        };
        /**
         * Notes the rows that have been selected within this pane and stores them internally
         *
         * @param notUpdating Whether the panes are updating themselves or not
         */
        SearchPane.prototype._updateSelection = function (notUpdating) {
            var _this = this;
            var settings = this.s.dt.settings()[0];
            var oApi = settings.oApi;
            var run = function () {
                _this.s.scrollTop = $$5(_this.s.dtPane.table().node()).parent()[0].scrollTop;
                if (_this.s.dt.page.info().serverSide && !_this.s.updating) {
                    if (!_this.s.serverSelecting) {
                        _this.s.serverSelect = _this.s.dtPane.rows({ selected: true }).data().toArray();
                        _this.s.dt.draw(false);
                    }
                }
                else if (notUpdating) {
                    _this._makeSelection();
                }
                oApi._fnProcessingDisplay(settings, false);
            };
            // If the processing display is enabled, we need to allow the browser
            // to draw it before performing our calculations
            if (settings.oFeatures.bProcessing) {
                oApi._fnProcessingDisplay(settings, true);
                setTimeout(run, 1);
            }
            else {
                run();
            }
        };
        /**
         * Takes in potentially undetected rows and adds them to the array if they are not yet featured
         *
         * @param filter the filter value of the potential row
         * @param display the display value of the potential row
         * @param sort the sort value of the potential row
         * @param type the type value of the potential row
         * @param arrayFilter the array to be populated
         * @param bins the bins to be populated
         */
        SearchPane.prototype._addOption = function (filter, display, sort, type, arrayFilter, bins) {
            // If the filter is an array then take a note of this, and add the elements to the arrayFilter array
            if (Array.isArray(filter) || filter instanceof dataTable$2.Api) {
                // Convert to an array so that we can work with it
                if (filter instanceof dataTable$2.Api) {
                    filter = filter.toArray();
                    display = display.toArray();
                }
                if (filter.length === display.length) {
                    for (var i = 0; i < filter.length; i++) {
                        // If we haven't seen this row before add it
                        if (!bins[filter[i]]) {
                            bins[filter[i]] = 1;
                            arrayFilter.push({
                                display: display[i],
                                filter: filter[i],
                                sort: sort[i],
                                type: type[i]
                            });
                        }
                        // Otherwise just increment the count
                        else {
                            bins[filter[i]]++;
                        }
                        this.s.rowData.totalOptions++;
                    }
                    return;
                }
                throw new Error('display and filter not the same length');
            }
            // If the values were affected by othogonal data and are not an array then check if it is already present
            else if (typeof this.s.colOpts.orthogonal === 'string') {
                if (!bins[filter]) {
                    bins[filter] = 1;
                    arrayFilter.push({
                        display: display,
                        filter: filter,
                        sort: sort,
                        type: type
                    });
                    this.s.rowData.totalOptions++;
                }
                else {
                    bins[filter]++;
                    this.s.rowData.totalOptions++;
                }
            }
            // Otherwise we must just be adding an option
            else {
                arrayFilter.push({
                    display: display,
                    filter: filter,
                    sort: sort,
                    type: type
                });
            }
        };
        /**
         * Method to construct the actual pane.
         *
         * @param selectedRows previously selected Rows to be reselected
         * @param dataIn Data that should be used to populate this pane
         * @param prevEl Reference to the previous element, used to ensure insert is in the correct location
         * @returns boolean to indicate whether this pane was the last one to have a selection made
         */
        SearchPane.prototype._buildPane = function (selectedRows, dataIn, prevEl) {
            var _this = this;
            if (selectedRows === void 0) { selectedRows = []; }
            if (dataIn === void 0) { dataIn = null; }
            if (prevEl === void 0) { prevEl = null; }
            // Aliases
            this.s.selections = [];
            // Other Variables
            var loadedFilter = this.s.dt.state.loaded();
            // If the listeners have not been set yet then using the latest state may result in funny errors
            if (this.s.listSet) {
                loadedFilter = this.s.dt.state();
            }
            // If it is not a custom pane in place
            if (this.s.colExists) {
                var idx = -1;
                if (loadedFilter && loadedFilter.searchPanes && loadedFilter.searchPanes.panes) {
                    for (var i = 0; i < loadedFilter.searchPanes.panes.length; i++) {
                        if (loadedFilter.searchPanes.panes[i].id === this.s.index) {
                            idx = i;
                            break;
                        }
                    }
                }
                // Perform checks that do not require populate pane to run
                if ((this.s.colOpts.show === false ||
                    this.s.colOpts.show !== undefined && this.s.colOpts.show !== true) &&
                    idx === -1) {
                    this.dom.container.addClass(this.classes.hidden);
                    this.s.displayed = false;
                    return false;
                }
                else if (this.s.colOpts.show === true || idx !== -1) {
                    this.s.displayed = true;
                }
                if (!this.s.dt.page.info().serverSide &&
                    (!dataIn ||
                        !dataIn.searchPanes ||
                        !dataIn.searchPanes.options)) {
                    // Only run populatePane if the data has not been collected yet
                    if (this.s.rowData.arrayFilter.length === 0) {
                        this.s.rowData.totalOptions = 0;
                        this._populatePane();
                        this.s.rowData.arrayOriginal = this.s.rowData.arrayFilter;
                        this.s.rowData.binsOriginal = this.s.rowData.bins;
                    }
                    var binLength = Object.keys(this.s.rowData.binsOriginal).length;
                    var uniqueRatio = this._uniqueRatio(binLength, this.s.dt.rows()[0].length);
                    // Don't show the pane if there isn't enough variance in the data, or there is only 1 entry
                    //  for that pane
                    if (this.s.displayed === false &&
                        ((this.s.colOpts.show === undefined && this.s.colOpts.threshold === null ?
                            uniqueRatio > this.c.threshold :
                            uniqueRatio > this.s.colOpts.threshold) ||
                            this.s.colOpts.show !== true && binLength <= 1)) {
                        this.dom.container.addClass(this.classes.hidden);
                        this.s.displayed = false;
                        return;
                    }
                    this.dom.container.addClass(this.classes.show);
                    this.s.displayed = true;
                }
                else if (dataIn && dataIn.searchPanes && dataIn.searchPanes.options) {
                    this._serverPopulate(dataIn);
                }
            }
            else {
                this.s.displayed = true;
            }
            // If the variance is accceptable then display the search pane
            this._displayPane();
            if (!this.s.listSet) {
                // Here, when the state is loaded if the data object on the original table is empty,
                //  then a state.clear() must have occurred, so delete all of the panes tables state objects too.
                this.dom.dtP.on('stateLoadParams.dtsp', function (e, settings, data) {
                    if ($$5.isEmptyObject(_this.s.dt.state.loaded())) {
                        $$5.each(data, function (index) {
                            delete data[index];
                        });
                    }
                });
            }
            // Add the container to the document in its original location
            if (prevEl !== null && this.dom.panesContainer.has(prevEl).length > 0) {
                this.dom.container.insertAfter(prevEl);
            }
            else {
                this.dom.panesContainer.prepend(this.dom.container);
            }
            // Declare the datatable for the pane
            var errMode = $$5.fn.dataTable.ext.errMode;
            $$5.fn.dataTable.ext.errMode = 'none';
            // eslint-disable-next-line no-extra-parens
            this.s.dtPane = this.dom.dtP.DataTable($$5.extend(true, this._getPaneConfig(), this.c.dtOpts, this.s.colOpts ? this.s.colOpts.dtOpts : {}, this.s.colOpts.options || !this.s.colExists ?
                {
                    createdRow: function (row, data) {
                        $$5(row).addClass(data.className);
                    }
                } :
                undefined, this.s.customPaneSettings !== null && this.s.customPaneSettings.dtOpts ?
                this.s.customPaneSettings.dtOpts :
                {}, $$5.fn.dataTable.versionCheck('2')
                ? {
                    layout: {
                        bottomLeft: null,
                        bottomRight: null,
                        topLeft: null,
                        topRight: null
                    }
                }
                : {}));
            this.dom.dtP.addClass(this.classes.table);
            // Getting column titles is a little messy
            var headerText = 'Custom Pane';
            if (this.s.customPaneSettings && this.s.customPaneSettings.header) {
                headerText = this.s.customPaneSettings.header;
            }
            else if (this.s.colOpts.header) {
                headerText = this.s.colOpts.header;
            }
            else if (this.s.colExists) {
                headerText = $$5.fn.dataTable.versionCheck('2')
                    ? this.s.dt.column(this.s.index).title()
                    : this.s.dt.settings()[0].aoColumns[this.s.index].sTitle;
            }
            headerText = this._escapeHTML(headerText);
            this.dom.searchBox.attr('placeholder', headerText);
            // As the pane table is not in the document yet we must initialise select ourselves
            // eslint-disable-next-line no-extra-parens
            $$5.fn.dataTable.select.init(this.s.dtPane);
            $$5.fn.dataTable.ext.errMode = errMode;
            // If it is not a custom pane
            if (this.s.colExists) {
                // Add all of the search options to the pane
                for (var i = 0, ien = this.s.rowData.arrayFilter.length; i < ien; i++) {
                    if (this.s.dt.page.info().serverSide) {
                        var row = this.addRow(this.s.rowData.arrayFilter[i].display, this.s.rowData.arrayFilter[i].filter, this.s.rowData.arrayFilter[i].sort, this.s.rowData.arrayFilter[i].type);
                        for (var _i = 0, _a = this.s.serverSelect; _i < _a.length; _i++) {
                            var option = _a[_i];
                            if (option.filter === this.s.rowData.arrayFilter[i].filter) {
                                this.s.serverSelecting = true;
                                row.select();
                                this.s.serverSelecting = false;
                            }
                        }
                    }
                    else if (!this.s.dt.page.info().serverSide && this.s.rowData.arrayFilter[i]) {
                        this.addRow(this.s.rowData.arrayFilter[i].display, this.s.rowData.arrayFilter[i].filter, this.s.rowData.arrayFilter[i].sort, this.s.rowData.arrayFilter[i].type);
                    }
                    else if (!this.s.dt.page.info().serverSide) {
                        // Just pass an empty string as the message will be calculated based on that in addRow()
                        this.addRow('', '', '', '');
                    }
                }
            }
            // eslint-disable-next-line no-extra-parens
            dataTable$2.select.init(this.s.dtPane);
            // If there are custom options set or it is a custom pane then get them
            if (this.s.colOpts.options ||
                this.s.customPaneSettings && this.s.customPaneSettings.options) {
                this._getComparisonRows();
            }
            // Display the pane
            this.s.dtPane.draw();
            this.s.dtPane.table().node().parentNode.scrollTop = this.s.scrollTop;
            this.adjustTopRow();
            this.setListeners();
            this.s.listSet = true;
            for (var _b = 0, selectedRows_1 = selectedRows; _b < selectedRows_1.length; _b++) {
                var selection = selectedRows_1[_b];
                if (selection) {
                    for (var _c = 0, _d = this.s.dtPane.rows().indexes().toArray(); _c < _d.length; _c++) {
                        var row = _d[_c];
                        if (this.s.dtPane.row(row).data() &&
                            selection.filter === this.s.dtPane.row(row).data().filter) {
                            // If this is happening when serverSide processing is happening then
                            //  different behaviour is needed
                            if (this.s.dt.page.info().serverSide) {
                                this.s.serverSelecting = true;
                                this.s.dtPane.row(row).select();
                                this.s.serverSelecting = false;
                            }
                            else {
                                this.s.dtPane.row(row).select();
                            }
                        }
                    }
                }
            }
            //  If SSP and the table is ready, apply the search for the pane
            if (this.s.dt.page.info().serverSide) {
                this.s.dtPane.search(this.dom.searchBox.val()).draw();
            }
            if ((this.c.initCollapsed && this.s.colOpts.initCollapsed !== false ||
                this.s.colOpts.initCollapsed) &&
                (this.c.collapse && this.s.colOpts.collapse !== false ||
                    this.s.colOpts.collapse)) {
                // If the pane has not initialised yet then we need to wait for it to do so before collapsing
                // Otherwise the container that the class is added to does not exist
                if (this.s.dtPane.settings()[0]._bInitComplete) {
                    this.collapse();
                }
                else {
                    this.s.dtPane.one('init', function () { return _this.collapse(); });
                }
            }
            // Reload the selection, searchbox entry and ordering from the previous state
            // Need to check here if SSP that this is the first draw, otherwise it will infinite loop
            if (loadedFilter &&
                loadedFilter.searchPanes &&
                loadedFilter.searchPanes.panes &&
                (!dataIn ||
                    dataIn.draw === 1)) {
                this._reloadSelect(loadedFilter);
                for (var _e = 0, _f = loadedFilter.searchPanes.panes; _e < _f.length; _e++) {
                    var pane = _f[_e];
                    if (pane.id === this.s.index) {
                        // Save some time by only triggering an input if there is a value
                        if (pane.searchTerm && pane.searchTerm.length > 0) {
                            this.dom.searchBox.val(pane.searchTerm).trigger('input');
                        }
                        if (pane.order) {
                            this.s.dtPane.order(pane.order).draw();
                        }
                        // Is the pane to be hidden or shown?
                        if (pane.collapsed) {
                            this.collapse();
                        }
                        else {
                            this.show();
                        }
                    }
                }
            }
            return true;
        };
        /**
         * Appends all of the HTML elements to their relevant parent Elements
         */
        SearchPane.prototype._displayPane = function () {
            // Empty everything to start again
            this.dom.dtP.empty();
            this.dom.topRow.empty().addClass(this.classes.topRow);
            // If there are more than 3 columns defined then make there be a smaller gap between the panes
            if (parseInt(this.c.layout.split('-')[1], 10) > 3) {
                this.dom.container.addClass(this.classes.smallGap);
            }
            this.dom.topRow
                .addClass(this.classes.subRowsContainer)
                .append(this.dom.upper.append(this.dom.searchCont))
                .append(this.dom.lower.append(this.dom.buttonGroup));
            // If no selections have been made in the pane then disable the clear button
            if (this.c.dtOpts.searching === false ||
                this.s.colOpts.dtOpts && this.s.colOpts.dtOpts.searching === false ||
                (!this.c.controls || !this.s.colOpts.controls) ||
                this.s.customPaneSettings &&
                    this.s.customPaneSettings.dtOpts &&
                    this.s.customPaneSettings.dtOpts.searching !== undefined &&
                    !this.s.customPaneSettings.dtOpts.searching) {
                this.dom.searchBox
                    .removeClass(this.classes.paneInputButton)
                    .addClass(this.classes.disabledButton)
                    .attr('disabled', 'true');
            }
            this.dom.searchBox.appendTo(this.dom.searchCont);
            // Create the contents of the searchCont div. Worth noting that this function will change when using semantic ui
            this._searchContSetup();
            // If the clear button is allowed to show then display it
            if (this.c.clear && this.c.controls && this.s.colOpts.controls) {
                this.dom.clear.appendTo(this.dom.buttonGroup);
            }
            if (this.c.orderable && this.s.colOpts.orderable && this.c.controls && this.s.colOpts.controls) {
                this.dom.nameButton.appendTo(this.dom.buttonGroup);
            }
            // If the count column is hidden then don't display the ordering button for it
            if (this.c.viewCount &&
                this.s.colOpts.viewCount &&
                this.c.orderable &&
                this.s.colOpts.orderable &&
                this.c.controls &&
                this.s.colOpts.controls) {
                this.dom.countButton.appendTo(this.dom.buttonGroup);
            }
            if ((this.c.collapse && this.s.colOpts.collapse !== false ||
                this.s.colOpts.collapse) &&
                this.c.controls && this.s.colOpts.controls) {
                this.dom.collapseButton.appendTo(this.dom.buttonGroup);
            }
            this.dom.container.prepend(this.dom.topRow).append(this.dom.dtP).show();
        };
        /**
         * Escape html characters within a string
         *
         * @param txt the string to be escaped
         * @returns the escaped string
         */
        SearchPane.prototype._escapeHTML = function (txt) {
            return txt
                .toString()
                .replace(/&amp;/g, '&')
                .replace(/&lt;/g, '<')
                .replace(/&gt;/g, '>')
                .replace(/&quot;/g, '"');
        };
        /**
         * Gets the options for the row for the customPanes
         *
         * @returns {object} The options for the row extended to include the options from the user.
         */
        SearchPane.prototype._getBonusOptions = function () {
            // We need to reset the thresholds as if they have a value in colOpts then that value will be used
            var defaultMutator = {
                threshold: null
            };
            return $$5.extend(true, {}, SearchPane.defaults, defaultMutator, this.c ? this.c : {});
        };
        /**
         * Gets the options for the row for the customPanes
         *
         * @returns {object} The options for the row extended to include the options from the user.
         */
        SearchPane.prototype._getOptions = function () {
            var table = this.s.dt;
            // We need to reset the thresholds as if they have a value in colOpts then that value will be used
            var defaultMutator = {
                collapse: null,
                emptyMessage: false,
                initCollapsed: null,
                threshold: null
            };
            var columnOptions = table.settings()[0].aoColumns[this.s.index].searchPanes;
            var colOpts = $$5.extend(true, {}, SearchPane.defaults, defaultMutator, columnOptions);
            if (columnOptions && columnOptions.hideCount && columnOptions.viewCount === undefined) {
                colOpts.viewCount = !columnOptions.hideCount;
            }
            return colOpts;
        };
        /**
         * Fill the array with the values that are currently being displayed in the table
         */
        SearchPane.prototype._populatePane = function () {
            this.s.rowData.arrayFilter = [];
            this.s.rowData.bins = {};
            var settings = this.s.dt.settings()[0];
            if (!this.s.dt.page.info().serverSide) {
                for (var _i = 0, _a = this.s.dt.rows().indexes().toArray(); _i < _a.length; _i++) {
                    var index = _a[_i];
                    this._populatePaneArray(index, this.s.rowData.arrayFilter, settings);
                }
            }
        };
        /**
         * This method decides whether a row should contribute to the pane or not
         *
         * @param filter the value that the row is to be filtered on
         * @param dataIndex the row index
         */
        SearchPane.prototype._search = function (filter, dataIndex) {
            var colOpts = this.s.colOpts;
            var table = this.s.dt;
            // For each item selected in the pane, check if it is available in the cell
            for (var _i = 0, _a = this.s.selections; _i < _a.length; _i++) {
                var colSelect = _a[_i];
                if (typeof colSelect === 'string' && typeof filter === 'string') {
                    // The filter value will not have the &amp; in place but a &,
                    // so we need to do a replace to make sure that they will match
                    colSelect = this._escapeHTML(colSelect);
                }
                // if the filter is an array then is the column present in it
                if (Array.isArray(filter)) {
                    if (colOpts.combiner === 'and') {
                        if (!filter.includes(colSelect)) {
                            return false;
                        }
                    }
                    else if (filter.includes(colSelect)) {
                        return true;
                    }
                }
                // if the filter is a function then does it meet the criteria of that function or not
                else if (typeof colSelect === 'function') {
                    if (colSelect.call(table, table.row(dataIndex).data(), dataIndex)) {
                        if (colOpts.combiner === 'or') {
                            return true;
                        }
                    }
                    // If the combiner is an "and" then we need to check against all possible selections
                    // so if it fails here then the and is not met and return false
                    else if (colOpts.combiner === 'and') {
                        return false;
                    }
                }
                // otherwise if the two filter values are equal then return true
                else if (filter === colSelect ||
                    // Loose type checking incase number type in column comparing to a string
                    // eslint-disable-next-line eqeqeq
                    !(typeof filter === 'string' && filter.length === 0) && filter == colSelect ||
                    colSelect === null && typeof filter === 'string' && filter === '') {
                    return true;
                }
            }
            // If the combiner is an and then we need to check against all possible selections
            // so return true here if so because it would have returned false earlier if it had failed
            if (colOpts.combiner === 'and') {
                return true;
            }
            // Otherwise it hasn't matched with anything by this point so it must be false
            return false;
        };
        /**
         * Creates the contents of the searchCont div
         *
         * NOTE This is overridden when semantic ui styling in order to integrate the search button into the text box.
         */
        SearchPane.prototype._searchContSetup = function () {
            if (this.c.controls && this.s.colOpts.controls) {
                this.dom.searchButton.appendTo(this.dom.searchLabelCont);
            }
            if (!(this.c.dtOpts.searching === false ||
                this.s.colOpts.dtOpts.searching === false ||
                this.s.customPaneSettings &&
                    this.s.customPaneSettings.dtOpts &&
                    this.s.customPaneSettings.dtOpts.searching !== undefined &&
                    !this.s.customPaneSettings.dtOpts.searching)) {
                this.dom.searchLabelCont.appendTo(this.dom.searchCont);
            }
        };
        /**
         * Adds outline to the pane when a selection has been made
         */
        SearchPane.prototype._searchExtras = function () {
            var updating = this.s.updating;
            this.s.updating = true;
            var filters = this.s.dtPane.rows({ selected: true }).data().pluck('filter').toArray();
            var nullIndex = filters.indexOf(this.emptyMessage());
            var container = $$5(this.s.dtPane.table().container());
            // If null index is found then search for empty cells as a filter.
            if (nullIndex > -1) {
                filters[nullIndex] = '';
            }
            // If a filter has been applied then outline the respective pane, remove it when it no longer is.
            if (filters.length > 0) {
                container.addClass(this.classes.selected);
            }
            else if (filters.length === 0) {
                container.removeClass(this.classes.selected);
            }
            this.s.updating = updating;
        };
        SearchPane.version = '2.1.2';
        SearchPane.classes = {
            bordered: 'dtsp-bordered',
            buttonGroup: 'dtsp-buttonGroup',
            buttonSub: 'dtsp-buttonSub',
            caret: 'dtsp-caret',
            clear: 'dtsp-clear',
            clearAll: 'dtsp-clearAll',
            clearButton: 'clearButton',
            collapseAll: 'dtsp-collapseAll',
            collapseButton: 'dtsp-collapseButton',
            container: 'dtsp-searchPane',
            countButton: 'dtsp-countButton',
            disabledButton: 'dtsp-disabledButton',
            hidden: 'dtsp-hidden',
            hide: 'dtsp-hide',
            layout: 'dtsp-',
            name: 'dtsp-name',
            nameButton: 'dtsp-nameButton',
            nameCont: 'dtsp-nameCont',
            narrow: 'dtsp-narrow',
            paneButton: 'dtsp-paneButton',
            paneInputButton: 'dtsp-paneInputButton',
            pill: 'dtsp-pill',
            rotated: 'dtsp-rotated',
            search: 'dtsp-search',
            searchCont: 'dtsp-searchCont',
            searchIcon: 'dtsp-searchIcon',
            searchLabelCont: 'dtsp-searchButtonCont',
            selected: 'dtsp-selected',
            smallGap: 'dtsp-smallGap',
            subRow1: 'dtsp-subRow1',
            subRow2: 'dtsp-subRow2',
            subRowsContainer: 'dtsp-subRowsContainer',
            title: 'dtsp-title',
            topRow: 'dtsp-topRow'
        };
        // Define SearchPanes default options
        SearchPane.defaults = {
            clear: true,
            collapse: true,
            combiner: 'or',
            container: function (dt) {
                return dt.table().container();
            },
            controls: true,
            dtOpts: {},
            emptyMessage: null,
            hideCount: false,
            i18n: {
                clearPane: '&times;',
                count: '{total}',
                emptyMessage: '<em>No data</em>'
            },
            initCollapsed: false,
            layout: 'auto',
            name: undefined,
            orderable: true,
            orthogonal: {
                display: 'display',
                filter: 'filter',
                hideCount: false,
                search: 'filter',
                show: undefined,
                sort: 'sort',
                threshold: 0.6,
                type: 'type',
                viewCount: true
            },
            preSelect: [],
            threshold: 0.6,
            viewCount: true
        };
        return SearchPane;
    }());

    var __extends$4 = (window && window.__extends) || (function () {
        var extendStatics = function (d, b) {
            extendStatics = Object.setPrototypeOf ||
                ({ __proto__: [] } instanceof Array && function (d, b) { d.__proto__ = b; }) ||
                function (d, b) { for (var p in b) if (b.hasOwnProperty(p)) d[p] = b[p]; };
            return extendStatics(d, b);
        };
        return function (d, b) {
            extendStatics(d, b);
            function __() { this.constructor = d; }
            d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
        };
    })();
    var SearchPaneST = /** @class */ (function (_super) {
        __extends$4(SearchPaneST, _super);
        function SearchPaneST(paneSettings, opts, index, panesContainer, panes) {
            return _super.call(this, paneSettings, opts, index, panesContainer, panes) || this;
        }
        /**
         * Populates the SearchPane based off of the data that has been recieved from the server
         *
         * This method overrides SearchPane's _serverPopulate() method
         *
         * @param dataIn The data that has been sent from the server
         */
        SearchPaneST.prototype._serverPopulate = function (dataIn) {
            this.s.rowData.binsShown = {};
            this.s.rowData.arrayFilter = [];
            if (dataIn.tableLength !== undefined) {
                this.s.tableLength = dataIn.tableLength;
                this.s.rowData.totalOptions = this.s.tableLength;
            }
            else if (this.s.tableLength === null || this.s.dt.rows()[0].length > this.s.tableLength) {
                this.s.tableLength = this.s.dt.rows()[0].length;
                this.s.rowData.totalOptions = this.s.tableLength;
            }
            var colTitle = this.s.dt.column(this.s.index).dataSrc();
            // If there is SP data for this column add it to the data array and bin
            if (dataIn.searchPanes.options[colTitle] !== undefined) {
                for (var _i = 0, _a = dataIn.searchPanes.options[colTitle]; _i < _a.length; _i++) {
                    var dataPoint = _a[_i];
                    this.s.rowData.arrayFilter.push({
                        display: dataPoint.label,
                        filter: dataPoint.value,
                        shown: +dataPoint.count,
                        sort: dataPoint.label,
                        total: +dataPoint.total,
                        type: dataPoint.label
                    });
                    this.s.rowData.binsShown[dataPoint.value] = +dataPoint.count;
                    this.s.rowData.bins[dataPoint.value] = +dataPoint.total;
                }
            }
            var binLength = Object.keys(this.s.rowData.bins).length;
            var uniqueRatio = this._uniqueRatio(binLength, this.s.tableLength);
            // Don't show the pane if there isnt enough variance in the data, or there is only 1 entry for that pane
            if (!this.s.colOpts.show &&
                this.s.displayed === false &&
                ((this.s.colOpts.show === undefined && this.s.colOpts.threshold === null ?
                    uniqueRatio > this.c.threshold :
                    uniqueRatio > this.s.colOpts.threshold) ||
                    this.s.colOpts.show !== true && binLength <= 1)) {
                this.dom.container.addClass(this.classes.hidden);
                this.s.displayed = false;
                return;
            }
            // Store the original data
            this.s.rowData.arrayOriginal = this.s.rowData.arrayFilter;
            this.s.rowData.binsOriginal = this.s.rowData.bins;
            // Flag this pane as being displayed
            this.s.displayed = true;
            // If the pane exists
            if (this.s.dtPane) {
                // Not the selections that have been made and remove all of the rows
                var selected = this.s.serverSelect;
                this.s.dtPane.rows().remove();
                // Add the rows that are to be shown into the pane
                for (var _b = 0, _c = this.s.rowData.arrayFilter; _b < _c.length; _b++) {
                    var data = _c[_b];
                    if (this._shouldAddRow(data)) {
                        var row = this.addRow(data.display, data.filter, data.sort, data.type);
                        // Select the row if it was selected before
                        for (var i = 0; i < selected.length; i++) {
                            var selection = selected[i];
                            if (selection.filter === data.filter) {
                                // This flag stops another request being made to the server
                                this.s.serverSelecting = true;
                                row.select();
                                this.s.serverSelecting = false;
                                // Remove the selection from the to select list and add it to the selected list
                                selected.splice(i, 1);
                                this.s.selections.push(data.filter);
                                break;
                            }
                        }
                    }
                }
                // Remake any selections that are no longer present
                for (var _d = 0, selected_1 = selected; _d < selected_1.length; _d++) {
                    var selection = selected_1[_d];
                    for (var _e = 0, _f = this.s.rowData.arrayOriginal; _e < _f.length; _e++) {
                        var data = _f[_e];
                        if (data.filter === selection.filter) {
                            var row = this.addRow(data.display, data.filter, data.sort, data.type);
                            this.s.serverSelecting = true;
                            row.select();
                            this.s.serverSelecting = false;
                            this.s.selections.push(data.filter);
                        }
                    }
                }
                // Store the selected rows
                this.s.serverSelect = this.s.dtPane.rows({ selected: true }).data().toArray();
                // Update the pane
                this.s.dtPane.draw();
            }
        };
        /**
         * This method updates the rows and their data within the SearchPanes
         *
         * SearchPaneCascade overrides this method
         */
        SearchPaneST.prototype.updateRows = function () {
            if (!this.s.dt.page.info().serverSide) {
                // Get the latest count values from the table
                this.s.rowData.binsShown = {};
                for (var _i = 0, _a = this.s.dt.rows({ search: 'applied' }).indexes().toArray(); _i < _a.length; _i++) {
                    var index = _a[_i];
                    this._updateShown(index, this.s.dt.settings()[0], this.s.rowData.binsShown);
                }
            }
            // Update the rows data to show the current counts
            for (var _b = 0, _c = this.s.dtPane.rows().data().toArray(); _b < _c.length; _b++) {
                var row = _c[_b];
                row.shown = typeof this.s.rowData.binsShown[row.filter] === 'number' ?
                    this.s.rowData.binsShown[row.filter] :
                    0;
                this.s.dtPane.row(row.index).data(row);
            }
            // Show updates in the pane
            this.s.dtPane.draw();
            this.s.dtPane.table().node().parentNode.scrollTop = this.s.scrollTop;
        };
        /**
         * Remove functionality from makeSelection - needs to be more advanced when tracking selections
         */
        SearchPaneST.prototype._makeSelection = function () {
            return;
        };
        /**
         * Blank method to remove reloading of selected rows - needs to be more advanced when tracking selections
         */
        SearchPaneST.prototype._reloadSelect = function () {
            return;
        };
        /**
         * Decides if a row should be added when being added from the server
         *
         * Overridden by SearchPaneCascade
         *
         * @param data the row data
         * @returns boolean indicating if the row should be added or not
         */
        // eslint-disable-next-line @typescript-eslint/no-unused-vars
        SearchPaneST.prototype._shouldAddRow = function (data) {
            return true;
        };
        /**
         * Updates the server selection list where appropriate
         */
        SearchPaneST.prototype._updateSelection = function () {
            if (this.s.dt.page.info().serverSide && !this.s.updating && !this.s.serverSelecting) {
                this.s.serverSelect = this.s.dtPane.rows({ selected: true }).data().toArray();
            }
        };
        /**
         * Used when binning the data for a column
         *
         * @param rowIdx The current row that is to be added to the bins
         * @param settings The datatables settings object
         * @param bins The bins object that is to be incremented
         */
        SearchPaneST.prototype._updateShown = function (rowIdx, settings, bins) {
            if (bins === void 0) { bins = this.s.rowData.binsShown; }
            var orth = typeof this.s.colOpts.orthogonal === 'string'
                ? this.s.colOpts.orthogonal
                : this.s.colOpts.orthogonal.search;
            var filter = settings.oApi._fnGetCellData(settings, rowIdx, this.s.index, orth);
            var add = function (f) {
                if (!bins[f]) {
                    bins[f] = 1;
                }
                else {
                    bins[f]++;
                }
            };
            if (Array.isArray(filter)) {
                for (var _i = 0, filter_1 = filter; _i < filter_1.length; _i++) {
                    var f = filter_1[_i];
                    add(f);
                }
            }
            else {
                add(filter);
            }
        };
        return SearchPaneST;
    }(SearchPane));

    var __extends$3 = (window && window.__extends) || (function () {
        var extendStatics = function (d, b) {
            extendStatics = Object.setPrototypeOf ||
                ({ __proto__: [] } instanceof Array && function (d, b) { d.__proto__ = b; }) ||
                function (d, b) { for (var p in b) if (b.hasOwnProperty(p)) d[p] = b[p]; };
            return extendStatics(d, b);
        };
        return function (d, b) {
            extendStatics(d, b);
            function __() { this.constructor = d; }
            d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
        };
    })();
    var $$4;
    function setJQuery$3(jq) {
        $$4 = jq;
    }
    var SearchPaneViewTotal = /** @class */ (function (_super) {
        __extends$3(SearchPaneViewTotal, _super);
        function SearchPaneViewTotal(paneSettings, opts, index, panesContainer, panes) {
            var _this = this;
            var override = {
                i18n: {
                    countFiltered: '{shown} ({total})'
                }
            };
            _this = _super.call(this, paneSettings, $$4.extend(override, opts), index, panesContainer, panes) || this;
            return _this;
        }
        /**
         * Gets the message that is to be used to indicate the count for each SearchPane row
         *
         * This method overrides _getMessage() in SearchPane and is overridden by SearchPaneCascadeViewTotal
         *
         * @param row The row object that is being processed
         * @returns string - the message that is to be shown for the count of each entry
         */
        SearchPaneViewTotal.prototype._getMessage = function (row) {
            var countMessage = this.s.dt.i18n('searchPanes.count', this.c.i18n.count);
            var filteredMessage = this.s.dt.i18n('searchPanes.countFiltered', this.c.i18n.countFiltered);
            return (this.s.filteringActive ? filteredMessage : countMessage)
                .replace(/{total}/g, row.total)
                .replace(/{shown}/g, row.shown);
        };
        /**
         * Overrides the blank method in SearchPane to return the number of times a given value is currently being displayed
         *
         * @param filter The filter value
         * @returns number - The number of times the value is shown
         */
        SearchPaneViewTotal.prototype._getShown = function (filter) {
            return this.s.rowData.binsShown && this.s.rowData.binsShown[filter] ?
                this.s.rowData.binsShown[filter] :
                0;
        };
        return SearchPaneViewTotal;
    }(SearchPaneST));

    var __extends$2 = (window && window.__extends) || (function () {
        var extendStatics = function (d, b) {
            extendStatics = Object.setPrototypeOf ||
                ({ __proto__: [] } instanceof Array && function (d, b) { d.__proto__ = b; }) ||
                function (d, b) { for (var p in b) if (b.hasOwnProperty(p)) d[p] = b[p]; };
            return extendStatics(d, b);
        };
        return function (d, b) {
            extendStatics(d, b);
            function __() { this.constructor = d; }
            d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
        };
    })();
    var $$3;
    function setJQuery$2(jq) {
        $$3 = jq;
    }
    var SearchPaneCascade = /** @class */ (function (_super) {
        __extends$2(SearchPaneCascade, _super);
        function SearchPaneCascade(paneSettings, opts, index, panesContainer, panes) {
            var _this = this;
            var override = {
                i18n: {
                    count: '{shown}'
                }
            };
            _this = _super.call(this, paneSettings, $$3.extend(override, opts), index, panesContainer, panes) || this;
            return _this;
        }
        /**
         * This method updates the rows and their data within the SearchPanes
         *
         * This overrides the method in SearchPane
         */
        SearchPaneCascade.prototype.updateRows = function () {
            // Note the currently selected values in the pane and remove all of the rows
            var selected = this.s.dtPane.rows({ selected: true }).data().toArray();
            if (this.s.colOpts.options ||
                this.s.customPaneSettings && this.s.customPaneSettings.options) {
                // If there are custom options set or it is a custom pane then get them
                this._getComparisonRows();
                var rows = this.s.dtPane.rows().toArray()[0];
                for (var i = 0; i < rows.length; i++) {
                    var row = this.s.dtPane.row(rows[i]);
                    var rowData = row.data();
                    if (rowData === undefined) {
                        continue;
                    }
                    if (rowData.shown === 0) {
                        row.remove();
                        rows = this.s.dtPane.rows().toArray()[0];
                        i--;
                        continue;
                    }
                    for (var _i = 0, selected_1 = selected; _i < selected_1.length; _i++) {
                        var selection = selected_1[_i];
                        if (rowData.filter === selection.filter) {
                            row.select();
                            selected.splice(i, 1);
                            this.s.selections.push(rowData.filter);
                            break;
                        }
                    }
                }
            }
            else {
                if (!this.s.dt.page.info().serverSide) {
                    // Get the latest count values from the table
                    this._activePopulatePane();
                    this.s.rowData.binsShown = {};
                    for (var _a = 0, _b = this.s.dt.rows({ search: 'applied' }).indexes().toArray(); _a < _b.length; _a++) {
                        var index = _b[_a];
                        this._updateShown(index, this.s.dt.settings()[0], this.s.rowData.binsShown);
                    }
                }
                this.s.dtPane.rows().remove();
                // Go over all of the rows that could be displayed
                for (var _c = 0, _d = this.s.rowData.arrayFilter; _c < _d.length; _c++) {
                    var data = _d[_c];
                    // Cascade - If there are no rows present in the table don't show the option
                    if (data.shown === 0) {
                        continue;
                    }
                    // Add the row to the pane
                    var row = this.addRow(data.display, data.filter, data.sort, data.type, undefined);
                    // Check if this row was selected
                    for (var i = 0; i < selected.length; i++) {
                        var selection = selected[i];
                        if (selection.filter === data.filter) {
                            row.select();
                            // Remove the row from the to find list and then add it to the selections list
                            selected.splice(i, 1);
                            this.s.selections.push(data.filter);
                            break;
                        }
                    }
                }
                // Add all of the rows that were previously selected but aren't any more
                for (var _e = 0, selected_2 = selected; _e < selected_2.length; _e++) {
                    var selection = selected_2[_e];
                    for (var _f = 0, _g = this.s.rowData.arrayOriginal; _f < _g.length; _f++) {
                        var data = _g[_f];
                        if (data.filter === selection.filter) {
                            var row = this.addRow(data.display, data.filter, data.sort, data.type, undefined);
                            row.select();
                            this.s.selections.push(data.filter);
                        }
                    }
                }
            }
            // Show updates in the pane
            this.s.dtPane.draw();
            this.s.dtPane.table().node().parentNode.scrollTop = this.s.scrollTop;
            // If client side updated the tables results
            if (!this.s.dt.page.info().serverSide) {
                this.s.dt.draw(false);
            }
        };
        /**
         * Fill the array with the values that are currently being displayed in the table
         */
        SearchPaneCascade.prototype._activePopulatePane = function () {
            this.s.rowData.arrayFilter = [];
            this.s.rowData.bins = {};
            var settings = this.s.dt.settings()[0];
            if (!this.s.dt.page.info().serverSide) {
                for (var _i = 0, _a = this.s.dt.rows({ search: 'applied' }).indexes().toArray(); _i < _a.length; _i++) {
                    var index = _a[_i];
                    this._populatePaneArray(index, this.s.rowData.arrayFilter, settings);
                }
            }
        };
        SearchPaneCascade.prototype._getComparisonRows = function () {
            // Find the appropriate options depending on whether this is a pane for a specific column or a custom pane
            var options = this.s.colOpts.options
                ? this.s.colOpts.options
                : this.s.customPaneSettings && this.s.customPaneSettings.options
                    ? this.s.customPaneSettings.options
                    : undefined;
            if (options === undefined) {
                return;
            }
            var allRows = this.s.dt.rows();
            var shownRows = this.s.dt.rows({ search: 'applied' });
            var tableValsTotal = allRows.data().toArray();
            var tableValsShown = shownRows.data().toArray();
            var rows = [];
            // Clear all of the other rows from the pane, only custom options are to be displayed when they are defined
            this.s.dtPane.clear();
            this.s.indexes = [];
            for (var _i = 0, options_1 = options; _i < options_1.length; _i++) {
                var comp = options_1[_i];
                // Initialise the object which is to be placed in the row
                var insert = comp.label !== '' ?
                    comp.label :
                    this.emptyMessage();
                var comparisonObj = {
                    className: comp.className,
                    display: insert,
                    filter: typeof comp.value === 'function' ? comp.value : [],
                    shown: 0,
                    sort: insert,
                    total: 0,
                    type: insert
                };
                // If a custom function is in place
                if (typeof comp.value === 'function') {
                    // Count the number of times the function evaluates to true for the original data in the Table
                    for (var i = 0; i < tableValsTotal.length; i++) {
                        if (comp.value.call(this.s.dt, tableValsTotal[i], allRows[0][i])) {
                            comparisonObj.total++;
                        }
                    }
                    for (var i = 0; i < tableValsShown.length; i++) {
                        if (comp.value.call(this.s.dt, tableValsShown[i], shownRows[0][i])) {
                            comparisonObj.shown++;
                        }
                    }
                    // Update the comparisonObj
                    if (typeof comparisonObj.filter !== 'function') {
                        comparisonObj.filter.push(comp.filter);
                    }
                }
                rows.push(this.addRow(comparisonObj.display, comparisonObj.filter, comparisonObj.sort, comparisonObj.type, comparisonObj.className, comparisonObj.total, comparisonObj.shown));
            }
            return rows;
        };
        /**
         * Gets the message that is to be used to indicate the count for each SearchPane row
         *
         * This method overrides _getMessage() in SearchPane and is overridden by SearchPaneCascadeViewTotal
         *
         * @param row The row object that is being processed
         * @returns string - the message that is to be shown for the count of each entry
         */
        SearchPaneCascade.prototype._getMessage = function (row) {
            return this.s.dt.i18n('searchPanes.count', this.c.i18n.count)
                .replace(/{total}/g, row.total)
                .replace(/{shown}/g, row.shown);
        };
        /**
         * Overrides the blank method in SearchPane to return the number of times a given value is currently being displayed
         *
         * @param filter The filter value
         * @returns number - The number of times the value is shown
         */
        SearchPaneCascade.prototype._getShown = function (filter) {
            return this.s.rowData.binsShown && this.s.rowData.binsShown[filter] ?
                this.s.rowData.binsShown[filter] :
                0;
        };
        /**
         * Decides if a row should be added when being added from the server
         *
         * Overrides method in by SearchPaneST
         *
         * @param data the row data
         * @returns boolean indicating if the row should be added or not
         */
        SearchPaneCascade.prototype._shouldAddRow = function (data) {
            return data.shown > 0;
        };
        return SearchPaneCascade;
    }(SearchPaneST));

    var __extends$1 = (window && window.__extends) || (function () {
        var extendStatics = function (d, b) {
            extendStatics = Object.setPrototypeOf ||
                ({ __proto__: [] } instanceof Array && function (d, b) { d.__proto__ = b; }) ||
                function (d, b) { for (var p in b) if (b.hasOwnProperty(p)) d[p] = b[p]; };
            return extendStatics(d, b);
        };
        return function (d, b) {
            extendStatics(d, b);
            function __() { this.constructor = d; }
            d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
        };
    })();
    var $$2;
    function setJQuery$1(jq) {
        $$2 = jq;
    }
    var SearchPaneCascadeViewTotal = /** @class */ (function (_super) {
        __extends$1(SearchPaneCascadeViewTotal, _super);
        function SearchPaneCascadeViewTotal(paneSettings, opts, index, panesContainer, panes) {
            var _this = this;
            var override = {
                i18n: {
                    count: '{total}',
                    countFiltered: '{shown} ({total})'
                }
            };
            _this = _super.call(this, paneSettings, $$2.extend(override, opts), index, panesContainer, panes) || this;
            return _this;
        }
        /**
         * Fill the array with the values that are currently being displayed in the table
         *
         * This method overrides _activePopulatePane() in SearchPaneCascade
         */
        SearchPaneCascadeViewTotal.prototype._activePopulatePane = function () {
            this.s.rowData.arrayFilter = [];
            this.s.rowData.binsShown = {};
            var settings = this.s.dt.settings()[0];
            if (!this.s.dt.page.info().serverSide) {
                for (var _i = 0, _a = this.s.dt.rows({ search: 'applied' }).indexes().toArray(); _i < _a.length; _i++) {
                    var index = _a[_i];
                    this._populatePaneArray(index, this.s.rowData.arrayFilter, settings, this.s.rowData.binsShown);
                }
            }
        };
        /**
         * Gets the message that is to be used to indicate the count for each SearchPane row
         *
         * This method overrides _getMessage() in SearchPaneCascade
         *
         * @param row The row object that is being processed
         * @returns string - the message that is to be shown for the count of each entry
         */
        SearchPaneCascadeViewTotal.prototype._getMessage = function (row) {
            var countMessage = this.s.dt.i18n('searchPanes.count', this.c.i18n.count);
            var filteredMessage = this.s.dt.i18n('searchPanes.countFiltered', this.c.i18n.countFiltered);
            return (this.s.filteringActive ? filteredMessage : countMessage)
                .replace(/{total}/g, row.total)
                .replace(/{shown}/g, row.shown);
        };
        return SearchPaneCascadeViewTotal;
    }(SearchPaneCascade));

    var $$1;
    var dataTable$1;
    function setJQuery(jq) {
        $$1 = jq;
        dataTable$1 = jq.fn.dataTable;
    }
    var SearchPanes = /** @class */ (function () {
        function SearchPanes(paneSettings, opts, fromPreInit, paneClass) {
            var _this = this;
            if (fromPreInit === void 0) { fromPreInit = false; }
            if (paneClass === void 0) { paneClass = SearchPane; }
            // Check that the required version of DataTables is included
            if (!dataTable$1 || !dataTable$1.versionCheck || !dataTable$1.versionCheck('1.10.0')) {
                throw new Error('SearchPane requires DataTables 1.10 or newer');
            }
            // Check that Select is included
            // eslint-disable-next-line no-extra-parens
            if (!dataTable$1.select) {
                throw new Error('SearchPane requires Select');
            }
            var table = new dataTable$1.Api(paneSettings);
            this.classes = $$1.extend(true, {}, SearchPanes.classes);
            // Get options from user
            this.c = $$1.extend(true, {}, SearchPanes.defaults, opts);
            // Add extra elements to DOM object including clear
            this.dom = {
                clearAll: $$1('<button type="button"/>')
                    .addClass(this.classes.clearAll)
                    .html(table.i18n('searchPanes.clearMessage', this.c.i18n.clearMessage)),
                collapseAll: $$1('<button type="button"/>')
                    .addClass(this.classes.collapseAll)
                    .html(table.i18n('searchPanes.collapseMessage', this.c.i18n.collapseMessage)),
                container: $$1('<div/>').addClass(this.classes.panes).html(table.i18n('searchPanes.loadMessage', this.c.i18n.loadMessage)),
                emptyMessage: $$1('<div/>').addClass(this.classes.emptyMessage),
                panes: $$1('<div/>').addClass(this.classes.container),
                showAll: $$1('<button type="button"/>')
                    .addClass(this.classes.showAll)
                    .addClass(this.classes.disabledButton)
                    .attr('disabled', 'true')
                    .html(table.i18n('searchPanes.showMessage', this.c.i18n.showMessage)),
                title: $$1('<div/>').addClass(this.classes.title),
                titleRow: $$1('<div/>').addClass(this.classes.titleRow)
            };
            this.s = {
                colOpts: [],
                dt: table,
                filterCount: 0,
                minPaneWidth: 260.0,
                page: 0,
                paging: false,
                pagingST: false,
                paneClass: paneClass,
                panes: [],
                selectionList: [],
                serverData: {},
                stateRead: false,
                updating: false
            };
            // Do not reinitialise if already initialised on table
            if (table.settings()[0]._searchPanes) {
                return;
            }
            this._getState();
            if (this.s.dt.page.info().serverSide) {
                var hostSettings = this.s.dt.settings()[0];
                // Listener to get the data into the server request before it is made
                this.s.dt.on('preXhr.dtsps', function (e, settings, data) {
                    if (hostSettings !== settings) {
                        return;
                    }
                    if (data.searchPanes === undefined) {
                        data.searchPanes = {};
                    }
                    if (data.searchPanes_null === undefined) {
                        data.searchPanes_null = {};
                    }
                    var src;
                    for (var _i = 0, _a = _this.s.selectionList; _i < _a.length; _i++) {
                        var selection = _a[_i];
                        src = _this.s.dt.column(selection.column).dataSrc();
                        if (data.searchPanes[src] === undefined) {
                            data.searchPanes[src] = {};
                        }
                        if (data.searchPanes_null[src] === undefined) {
                            data.searchPanes_null[src] = {};
                        }
                        for (var i = 0; i < selection.rows.length; i++) {
                            data.searchPanes[src][i] = selection.rows[i];
                            if (data.searchPanes[src][i] === null) {
                                data.searchPanes_null[src][i] = true;
                            }
                        }
                    }
                    if (_this.s.selectionList.length > 0) {
                        data.searchPanesLast = src;
                    }
                    // Config options that will change how the querying is done
                    data.searchPanes_options = {
                        cascade: _this.c.cascadePanes,
                        viewCount: _this.c.viewCount,
                        viewTotal: _this.c.viewTotal
                    };
                });
            }
            this._setXHR();
            table.settings()[0]._searchPanes = this;
            if (this.s.dt.settings()[0]._bInitComplete || fromPreInit) {
                this._paneDeclare(table, paneSettings, opts);
            }
            else {
                table.one('preInit.dtsps', function () {
                    _this._paneDeclare(table, paneSettings, opts);
                });
            }
            return this;
        }
        /**
         * Clear the selections of all of the panes
         */
        SearchPanes.prototype.clearSelections = function () {
            for (var _i = 0, _a = this.s.panes; _i < _a.length; _i++) {
                var pane = _a[_i];
                if (pane.s.dtPane) {
                    pane.s.scrollTop = pane.s.dtPane.table().node().parentNode.scrollTop;
                }
            }
            // Load in all of the searchBoxes in the documents
            var searches = this.dom.container.find('.' + this.classes.search.replace(/\s+/g, '.'));
            // For each searchBox set the input text to be empty and then trigger
            // an input on them so that they no longer filter the panes
            searches.each(function () {
                $$1(this).val('').trigger('input');
            });
            // Clear the selectionList
            this.s.selectionList = [];
            var returnArray = [];
            for (var _b = 0, _c = this.s.panes; _b < _c.length; _b++) {
                var pane = _c[_b];
                if (pane.s.dtPane) {
                    returnArray.push(pane.clearPane());
                }
            }
            return returnArray;
        };
        /**
         * returns the container node for the searchPanes
         */
        SearchPanes.prototype.getNode = function () {
            return this.dom.container;
        };
        /**
         * rebuilds all of the panes
         */
        SearchPanes.prototype.rebuild = function (targetIdx, maintainSelection) {
            if (targetIdx === void 0) { targetIdx = false; }
            if (maintainSelection === void 0) { maintainSelection = false; }
            this.dom.emptyMessage.detach();
            // As a rebuild from scratch is required, empty the searchpanes container.
            if (targetIdx === false) {
                this.dom.panes.empty();
            }
            // Rebuild each pane individually, if a specific pane has been selected then only rebuild that one
            var returnArray = [];
            for (var _i = 0, _a = this.s.panes; _i < _a.length; _i++) {
                var pane = _a[_i];
                if (targetIdx === false || pane.s.index === targetIdx) {
                    pane.clearData();
                    pane.rebuildPane(this.s.dt.page.info().serverSide ?
                        this.s.serverData :
                        undefined, maintainSelection);
                    this.dom.panes.append(pane.dom.container);
                    returnArray.push(pane);
                }
            }
            this._updateSelection();
            // Attach panes, clear buttons, and title bar to the document
            this._updateFilterCount();
            this._attachPaneContainer();
            this._initSelectionListeners(false);
            // If the selections are to be maintained, then it is safe to assume that paging is also to be maintained
            // Otherwise, the paging should be reset
            this.s.dt.draw(!maintainSelection);
            // Resize the panes incase there has been a change
            this.resizePanes();
            // If a single pane has been rebuilt then return only that pane
            return returnArray.length === 1 ? returnArray[0] : returnArray;
        };
        /**
         * Resizes all of the panes
         */
        SearchPanes.prototype.resizePanes = function () {
            if (this.c.layout === 'auto') {
                var contWidth = $$1(this.s.dt.searchPanes.container()).width();
                var target = Math.floor(contWidth / this.s.minPaneWidth); // The neatest number of panes per row
                var highest_1 = 1;
                var highestmod_1 = 0;
                // Get the indexes of all of the displayed panes
                var dispIndex = [];
                for (var _i = 0, _a = this.s.panes; _i < _a.length; _i++) {
                    var pane = _a[_i];
                    if (pane.s.displayed) {
                        dispIndex.push(pane.s.index);
                    }
                }
                var displayCount = dispIndex.length;
                // If the neatest number is the number we have then use this.
                if (target === displayCount) {
                    highest_1 = target;
                }
                else {
                    // Go from the target down and find the value with the most panes left over, this will be the best fit
                    for (var ppr = target; ppr > 1; ppr--) {
                        var rem = displayCount % ppr;
                        if (rem === 0) {
                            highest_1 = ppr;
                            highestmod_1 = 0;
                            break;
                        }
                        // If there are more left over at this amount of panes per row (ppr)
                        // then it fits better so new values
                        else if (rem > highestmod_1) {
                            highest_1 = ppr;
                            highestmod_1 = rem;
                        }
                    }
                }
                // If there is a perfect fit then none are to be wider
                var widerIndexes_1 = highestmod_1 !== 0 ? dispIndex.slice(dispIndex.length - highestmod_1, dispIndex.length) : [];
                this.s.panes.forEach(function (pane) {
                    // Resize the pane with the new layout
                    if (pane.s.displayed) {
                        pane.resize('columns-' + (!widerIndexes_1.includes(pane.s.index) ? highest_1 : highestmod_1));
                    }
                });
            }
            else {
                for (var _b = 0, _c = this.s.panes; _b < _c.length; _b++) {
                    var pane = _c[_b];
                    pane.adjustTopRow();
                }
            }
            return this;
        };
        /**
         * Holder method that is userd in SearchPanesST to set listeners that have an effect on other panes
         *
         * @param isPreselect boolean to indicate if the preselect array is to override the current selection list.
         */
        SearchPanes.prototype._initSelectionListeners = function (isPreselect) {
            return;
        };
        /**
         * Blank method that is overridden in SearchPanesST to retrieve the totals from the server data
         */
        SearchPanes.prototype._serverTotals = function () {
            return;
        };
        /**
         * Set's the xhr listener so that SP can extact appropriate data from the response
         */
        SearchPanes.prototype._setXHR = function () {
            var _this = this;
            var hostSettings = this.s.dt.settings()[0];
            var run = function (json) {
                if (json && json.searchPanes && json.searchPanes.options) {
                    _this.s.serverData = json;
                    _this.s.serverData.tableLength = json.recordsTotal;
                    _this._serverTotals();
                }
            };
            // We are using the xhr event to rebuild the panes if required due to viewTotal being enabled
            // If viewTotal is not enabled then we simply update the data from the server
            this.s.dt.on('xhr.dtsps', function (e, settings, json) {
                if (hostSettings === settings) {
                    run(json);
                }
            });
            // Account for the initial JSON fetch having already completed
            run(this.s.dt.ajax.json());
        };
        /**
         * Set's the function that is to be performed when a state is loaded
         *
         * Overridden by the method in SearchPanesST
         */
        SearchPanes.prototype._stateLoadListener = function () {
            var _this = this;
            var hostSettings = this.s.dt.settings()[0];
            this.s.dt.on('stateLoadParams.dtsps', function (e, settings, data) {
                if (data.searchPanes === undefined || settings !== hostSettings) {
                    return;
                }
                _this.clearSelections();
                // Set the selection list for the panes so that the correct
                // rows can be reselected and in the right order
                _this.s.selectionList =
                    data.searchPanes.selectionList ?
                        data.searchPanes.selectionList :
                        [];
                // Find the panes that match from the state and the actual instance
                if (data.searchPanes.panes) {
                    for (var _i = 0, _a = data.searchPanes.panes; _i < _a.length; _i++) {
                        var loadedPane = _a[_i];
                        for (var _b = 0, _c = _this.s.panes; _b < _c.length; _b++) {
                            var pane = _c[_b];
                            if (loadedPane.id === pane.s.index && pane.s.dtPane) {
                                // Set the value of the searchbox
                                pane.dom.searchBox.val(loadedPane.searchTerm);
                                // Set the value of the order
                                pane.s.dtPane.order(loadedPane.order);
                            }
                        }
                    }
                }
                _this._makeSelections(_this.s.selectionList);
            });
        };
        /**
         * Updates the selectionList when cascade is not in place
         *
         * Overridden in SearchPanesST
         */
        SearchPanes.prototype._updateSelection = function () {
            this.s.selectionList = [];
            for (var _i = 0, _a = this.s.panes; _i < _a.length; _i++) {
                var pane = _a[_i];
                if (pane.s.dtPane) {
                    var rows = pane.s.dtPane.rows({ selected: true }).data().toArray().map(function (el) { return el.filter; });
                    if (rows.length) {
                        this.s.selectionList.push({
                            column: pane.s.index,
                            rows: rows
                        });
                    }
                }
            }
        };
        /**
         * Attach the panes, buttons and title to the document
         */
        SearchPanes.prototype._attach = function () {
            var _this = this;
            this.dom.titleRow
                .removeClass(this.classes.hide)
                .detach()
                .append(this.dom.title);
            // If the clear button is permitted attach it
            if (this.c.clear) {
                this.dom.clearAll
                    .appendTo(this.dom.titleRow)
                    .on('click.dtsps', function () { return _this.clearSelections(); });
            }
            if (this.c.collapse) {
                this.dom.showAll.appendTo(this.dom.titleRow);
                this.dom.collapseAll.appendTo(this.dom.titleRow);
                this._setCollapseListener();
            }
            // Attach the container for each individual pane to the overall container
            for (var _i = 0, _a = this.s.panes; _i < _a.length; _i++) {
                var pane = _a[_i];
                this.dom.panes.append(pane.dom.container);
            }
            // Attach everything to the document
            this.dom.container
                .text('')
                .removeClass(this.classes.hide)
                .append(this.dom.titleRow)
                .append(this.dom.panes);
            // WORKAROUND
            this.s.panes.forEach(function (pane) { return pane.setListeners(); });
            if ($$1('div.' + this.classes.container).length === 0) {
                this.dom.container.prependTo(this.s.dt);
            }
        };
        /**
         * If there are no panes to display then this method is called to either
         * display a message in their place or hide them completely.
         */
        SearchPanes.prototype._attachMessage = function () {
            // Create a message to display on the screen
            var message;
            try {
                message = this.s.dt.i18n('searchPanes.emptyPanes', this.c.i18n.emptyPanes);
            }
            catch (error) {
                message = null;
            }
            // If the message is an empty string then searchPanes.emptyPanes is undefined,
            // therefore the pane container should be removed from the display
            if (message === null) {
                this.dom.container.addClass(this.classes.hide);
                this.dom.titleRow.removeClass(this.classes.hide);
                return;
            }
            // Otherwise display the message
            this.dom.container.removeClass(this.classes.hide);
            this.dom.titleRow.addClass(this.classes.hide);
            this.dom.emptyMessage.html(message).appendTo(this.dom.container);
        };
        /**
         * Attaches the panes to the document and displays a message or hides if there are none
         */
        SearchPanes.prototype._attachPaneContainer = function () {
            // If a pane is to be displayed then attach the normal pane output
            for (var _i = 0, _a = this.s.panes; _i < _a.length; _i++) {
                var pane = _a[_i];
                if (pane.s.displayed === true) {
                    this._attach();
                    return;
                }
            }
            // Otherwise attach the custom message or remove the container from the display
            this._attachMessage();
        };
        /**
         * Checks which panes are collapsed and then performs relevant actions to the collapse/show all buttons
         */
        SearchPanes.prototype._checkCollapse = function () {
            var disableClose = true;
            var disableShow = true;
            for (var _i = 0, _a = this.s.panes; _i < _a.length; _i++) {
                var pane = _a[_i];
                if (pane.s.displayed) {
                    // If the pane is not collapsed
                    if (!pane.dom.collapseButton.hasClass(pane.classes.rotated)) {
                        // Enable the collapse all button
                        this.dom.collapseAll.removeClass(this.classes.disabledButton).removeAttr('disabled');
                        disableClose = false;
                    }
                    else {
                        // Otherwise enable the show all button
                        this.dom.showAll.removeClass(this.classes.disabledButton).removeAttr('disabled');
                        disableShow = false;
                    }
                }
            }
            // If this flag is still true, no panes are open so the close button should be disabled
            if (disableClose) {
                this.dom.collapseAll.addClass(this.classes.disabledButton).attr('disabled', 'true');
            }
            // If this flag is still true, no panes are closed so the show button should be disabled
            if (disableShow) {
                this.dom.showAll.addClass(this.classes.disabledButton).attr('disabled', 'true');
            }
        };
        /**
         * Attaches the message to the document but does not add any panes
         */
        SearchPanes.prototype._checkMessage = function () {
            // If a pane is to be displayed then attach the normal pane output
            for (var _i = 0, _a = this.s.panes; _i < _a.length; _i++) {
                var pane = _a[_i];
                if (pane.s.displayed === true) {
                    // Ensure that the empty message is removed if a pane is displayed
                    this.dom.emptyMessage.detach();
                    this.dom.titleRow.removeClass(this.classes.hide);
                    return;
                }
            }
            // Otherwise attach the custom message or remove the container from the display
            this._attachMessage();
        };
        /**
         * Collapses all of the panes
         */
        SearchPanes.prototype._collapseAll = function () {
            for (var _i = 0, _a = this.s.panes; _i < _a.length; _i++) {
                var pane = _a[_i];
                pane.collapse();
            }
        };
        /**
         * Finds a pane based upon the name of that pane
         *
         * @param name string representing the name of the pane
         * @returns SearchPane The pane which has that name
         */
        SearchPanes.prototype._findPane = function (name) {
            for (var _i = 0, _a = this.s.panes; _i < _a.length; _i++) {
                var pane = _a[_i];
                if (name === pane.s.name) {
                    return pane;
                }
            }
        };
        /**
         * Gets the selection list from the previous state and stores it in the selectionList Property
         */
        SearchPanes.prototype._getState = function () {
            var loadedFilter = this.s.dt.state.loaded();
            if (loadedFilter && loadedFilter.searchPanes && loadedFilter.searchPanes.selectionList) {
                this.s.selectionList = loadedFilter.searchPanes.selectionList;
            }
        };
        SearchPanes.prototype._makeSelections = function (selectList) {
            for (var _i = 0, selectList_1 = selectList; _i < selectList_1.length; _i++) {
                var selection = selectList_1[_i];
                var pane = void 0;
                for (var _a = 0, _b = this.s.panes; _a < _b.length; _a++) {
                    var p = _b[_a];
                    if (p.s.index === selection.column) {
                        pane = p;
                        break;
                    }
                }
                if (pane && pane.s.dtPane) {
                    for (var j = 0; j < pane.s.dtPane.rows().data().toArray().length; j++) {
                        if (selection.rows.includes(typeof pane.s.dtPane.row(j).data().filter === 'function' ?
                            pane.s.dtPane.cell(j, 0).data() :
                            pane.s.dtPane.row(j).data().filter)) {
                            pane.s.dtPane.row(j).select();
                        }
                    }
                    pane.updateTable();
                }
            }
        };
        /**
         * Declares the instances of individual searchpanes dependant on the number of columns.
         * It is necessary to run this once preInit has completed otherwise no panes will be
         * created as the column count will be 0.
         *
         * @param table the DataTable api for the parent table
         * @param paneSettings the settings passed into the constructor
         * @param opts the options passed into the constructor
         */
        SearchPanes.prototype._paneDeclare = function (table, paneSettings, opts) {
            var _this = this;
            // Create Panes
            table
                .columns(this.c.columns.length > 0 ? this.c.columns : undefined)
                .eq(0)
                .each(function (idx) {
                _this.s.panes.push(new _this.s.paneClass(paneSettings, opts, idx, _this.dom.panes));
            });
            // If there is any extra custom panes defined then create panes for them too
            var colCount = table.columns().eq(0).toArray().length;
            for (var i = 0; i < this.c.panes.length; i++) {
                var id = colCount + i;
                this.s.panes.push(new this.s.paneClass(paneSettings, opts, id, this.dom.panes, this.c.panes[i]));
            }
            // If a custom ordering is being used
            if (this.c.order.length > 0) {
                // Make a new Array of panes based upon the order
                this.s.panes = this.c.order.map(function (name) { return _this._findPane(name); });
            }
            // If this internal property is true then the DataTable has been initialised already
            if (this.s.dt.settings()[0]._bInitComplete) {
                this._startup(table);
            }
            else {
                // Otherwise add the paneStartup function to the list of functions
                // that are to be run when the table is initialised. This will garauntee that the
                // panes are initialised before the init event and init Complete callback is fired
                this.s.dt.settings()[0].aoInitComplete.push({
                    fn: function () { return _this._startup(table); }
                });
            }
        };
        /**
         * Sets the listeners for the collapse and show all buttons
         * Also sets and performs checks on current panes to see if they are collapsed
         */
        SearchPanes.prototype._setCollapseListener = function () {
            var _this = this;
            this.dom.collapseAll.on('click.dtsps', function () {
                _this._collapseAll();
                _this.dom.collapseAll.addClass(_this.classes.disabledButton).attr('disabled', 'true');
                _this.dom.showAll.removeClass(_this.classes.disabledButton).removeAttr('disabled');
                _this.s.dt.state.save();
            });
            this.dom.showAll.on('click.dtsps', function () {
                _this._showAll();
                _this.dom.showAll.addClass(_this.classes.disabledButton).attr('disabled', 'true');
                _this.dom.collapseAll.removeClass(_this.classes.disabledButton).removeAttr('disabled');
                _this.s.dt.state.save();
            });
            for (var _i = 0, _a = this.s.panes; _i < _a.length; _i++) {
                var pane = _a[_i];
                // We want to make the same check whenever there is a collapse/expand
                pane.dom.topRow.on('collapse.dtsps', function () { return _this._checkCollapse(); });
            }
            this._checkCollapse();
        };
        /**
         * Shows all of the panes
         */
        SearchPanes.prototype._showAll = function () {
            for (var _i = 0, _a = this.s.panes; _i < _a.length; _i++) {
                var pane = _a[_i];
                pane.show();
            }
        };
        /**
         * Initialises the tables previous/preset selections and initialises callbacks for events
         *
         * @param table the parent table for which the searchPanes are being created
         */
        SearchPanes.prototype._startup = function (table) {
            var _this = this;
            // Attach clear button and title bar to the document
            this._attach();
            this.dom.panes.empty();
            var hostSettings = this.s.dt.settings()[0];
            for (var _i = 0, _a = this.s.panes; _i < _a.length; _i++) {
                var pane = _a[_i];
                pane.rebuildPane(Object.keys(this.s.serverData).length > 0 ? this.s.serverData : undefined);
                this.dom.panes.append(pane.dom.container);
            }
            // If the layout is set to auto then the panes need to be resized to their best fit
            if (this.c.layout === 'auto') {
                this.resizePanes();
            }
            var loadedFilter = this.s.dt.state.loaded();
            // Reset the paging if that has been saved in the state
            if (!this.s.stateRead && loadedFilter) {
                this.s.dt
                    .page(loadedFilter.start / this.s.dt.page.len())
                    .draw('page');
            }
            this.s.stateRead = true;
            this._checkMessage();
            // When a draw is called on the DataTable, update all of the panes incase the data in the DataTable has changed
            table.on('preDraw.dtsps', function () {
                // Check that the panes are not updating to avoid infinite loops
                // Also check that this draw is not due to paging
                if (!_this.s.updating && !_this.s.paging) {
                    _this._updateFilterCount();
                    _this._updateSelection();
                }
                // Paging flag reset - we only need to dodge the draw once
                _this.s.paging = false;
            });
            $$1(window).on('resize.dtsps', dataTable$1.util.throttle(function () { return _this.resizePanes(); }));
            // Whenever a state save occurs store the selection list in the state object
            this.s.dt.on('stateSaveParams.dtsps', function (e, settings, data) {
                if (settings !== hostSettings) {
                    return;
                }
                if (data.searchPanes === undefined) {
                    data.searchPanes = {};
                }
                data.searchPanes.selectionList = _this.s.selectionList;
            });
            this._stateLoadListener();
            // Listener for paging on main table
            table.off('page.dtsps page-nc.dtsps').on('page.dtsps page-nc.dtsps', function (e, s) {
                _this.s.paging = true;
                // This is an indicator to any selection tracking classes that paging has occured
                // It has to happen here so that we don't stack event listeners unnecessarily
                // The value is only ever set back to false in the SearchPanesST class
                // Equally it is never read in this class
                _this.s.pagingST = true;
                _this.s.page = _this.s.dt.page();
            });
            if (this.s.dt.page.info().serverSide) {
                table.off('preXhr.dtsps').on('preXhr.dtsps', function (e, settings, data) {
                    if (settings !== hostSettings) {
                        return;
                    }
                    if (!data.searchPanes) {
                        data.searchPanes = {};
                    }
                    if (!data.searchPanes_null) {
                        data.searchPanes_null = {};
                    }
                    // Count how many filters are being applied
                    var filterCount = 0;
                    for (var _i = 0, _a = _this.s.panes; _i < _a.length; _i++) {
                        var pane = _a[_i];
                        var src = _this.s.dt.column(pane.s.index).dataSrc();
                        if (!data.searchPanes[src]) {
                            data.searchPanes[src] = {};
                        }
                        if (!data.searchPanes_null[src]) {
                            data.searchPanes_null[src] = {};
                        }
                        if (pane.s.dtPane) {
                            var rowData = pane.s.dtPane.rows({ selected: true }).data().toArray();
                            for (var i = 0; i < rowData.length; i++) {
                                data.searchPanes[src][i] = rowData[i].filter;
                                if (!data.searchPanes[src][i]) {
                                    data.searchPanes_null[src][i] = true;
                                }
                                filterCount++;
                            }
                        }
                    }
                    // If there is a filter to be applied, then we need to read from the start of the result set
                    // and set the paging to 0. This matches the behaviour of client side processing
                    if (filterCount > 0) {
                        // If the number of filters has changed we need to read from the start of the
                        // result set and reset the paging
                        if (filterCount !== _this.s.filterCount) {
                            data.start = 0;
                            _this.s.page = 0;
                        }
                        // Otherwise it is a paging request and we need to read from whatever the paging has been set to
                        else {
                            data.start = _this.s.page * _this.s.dt.page.len();
                        }
                        _this.s.dt.page(_this.s.page);
                        _this.s.filterCount = filterCount;
                    }
                    if (_this.s.selectionList.length > 0) {
                        data.searchPanesLast = _this.s.dt
                            .column(_this.s.selectionList[_this.s.selectionList.length - 1].column)
                            .dataSrc();
                    }
                    // Config options that will change how the querying is done
                    data.searchPanes_options = {
                        cascade: _this.c.cascadePanes,
                        viewCount: _this.c.viewCount,
                        viewTotal: _this.c.viewTotal
                    };
                });
            }
            else {
                table.on('preXhr.dtsps', function () { return _this.s.panes.forEach(function (pane) { return pane.clearData(); }); });
            }
            // If the data is reloaded from the server then it is possible that it has changed completely,
            // so we need to rebuild the panes
            this.s.dt.on('xhr.dtsps', function (e, settings) {
                if (settings.nTable !== _this.s.dt.table().node()) {
                    return;
                }
                if (!_this.s.dt.page.info().serverSide) {
                    var processing_1 = false;
                    _this.s.dt.one('preDraw.dtsps', function () {
                        if (processing_1) {
                            return;
                        }
                        var page = _this.s.dt.page();
                        processing_1 = true;
                        _this.s.updating = true;
                        _this.dom.panes.empty();
                        for (var _i = 0, _a = _this.s.panes; _i < _a.length; _i++) {
                            var pane = _a[_i];
                            pane.clearData(); // Clears all of the bins and will mean that the data has to be re-read
                            // Pass a boolean to say whether this is the last choice made for maintaining selections
                            // when rebuilding
                            pane.rebuildPane(undefined, true);
                            _this.dom.panes.append(pane.dom.container);
                        }
                        if (!_this.s.dt.page.info().serverSide) {
                            _this.s.dt.draw();
                        }
                        _this.s.updating = false;
                        _this._updateSelection();
                        _this._checkMessage();
                        _this.s.dt.one('draw.dtsps', function () {
                            _this.s.updating = true;
                            _this.s.dt.page(page).draw(false);
                            _this.s.updating = false;
                        });
                    });
                }
            });
            // PreSelect any selections which have been defined using the preSelect option
            var selectList = this.c.preSelect;
            if (loadedFilter && loadedFilter.searchPanes && loadedFilter.searchPanes.selectionList) {
                selectList = loadedFilter.searchPanes.selectionList;
            }
            this._makeSelections(selectList);
            // Update the title bar to show how many filters have been selected
            this._updateFilterCount();
            // If the table is destroyed and restarted then clear the selections so that they do not persist.
            table.on('destroy.dtsps', function (e, settings) {
                if (settings !== hostSettings) {
                    return;
                }
                for (var _i = 0, _a = _this.s.panes; _i < _a.length; _i++) {
                    var pane = _a[_i];
                    pane.destroy();
                }
                table.off('.dtsps');
                _this.dom.showAll.off('.dtsps');
                _this.dom.clearAll.off('.dtsps');
                _this.dom.collapseAll.off('.dtsps');
                $$1(table.table().node()).off('.dtsps');
                _this.dom.container.detach();
                _this.clearSelections();
            });
            if (this.c.collapse) {
                this._setCollapseListener();
            }
            // When the clear All button has been pressed clear all of the selections in the panes
            if (this.c.clear) {
                this.dom.clearAll.on('click.dtsps', function () { return _this.clearSelections(); });
            }
            hostSettings._searchPanes = this;
            // This state save is required so that state is maintained over multiple refreshes if no actions are made
            this.s.dt.state.save();
        };
        /**
         * Updates the number of filters that have been applied in the title
         */
        SearchPanes.prototype._updateFilterCount = function () {
            var filterCount = 0;
            // Add the number of all of the filters throughout the panes
            for (var _i = 0, _a = this.s.panes; _i < _a.length; _i++) {
                var pane = _a[_i];
                if (pane.s.dtPane) {
                    filterCount += pane.getPaneCount();
                }
            }
            // Run the message through the internationalisation method to improve readability
            this.dom.title.html(this.s.dt.i18n('searchPanes.title', this.c.i18n.title, filterCount));
            if (this.c.filterChanged && typeof this.c.filterChanged === 'function') {
                this.c.filterChanged.call(this.s.dt, filterCount);
            }
            if (filterCount === 0) {
                this.dom.clearAll.addClass(this.classes.disabledButton).attr('disabled', 'true');
            }
            else {
                this.dom.clearAll.removeClass(this.classes.disabledButton).removeAttr('disabled');
            }
        };
        SearchPanes.version = '2.1.2';
        SearchPanes.classes = {
            clear: 'dtsp-clear',
            clearAll: 'dtsp-clearAll',
            collapseAll: 'dtsp-collapseAll',
            container: 'dtsp-searchPanes',
            disabledButton: 'dtsp-disabledButton',
            emptyMessage: 'dtsp-emptyMessage',
            hide: 'dtsp-hidden',
            panes: 'dtsp-panesContainer',
            search: 'dtsp-search',
            showAll: 'dtsp-showAll',
            title: 'dtsp-title',
            titleRow: 'dtsp-titleRow'
        };
        // Define SearchPanes default options
        SearchPanes.defaults = {
            cascadePanes: false,
            clear: true,
            collapse: true,
            columns: [],
            container: function (dt) {
                return dt.table().container();
            },
            filterChanged: undefined,
            i18n: {
                clearMessage: 'Clear All',
                clearPane: '&times;',
                collapse: {
                    0: 'SearchPanes',
                    _: 'SearchPanes (%d)'
                },
                collapseMessage: 'Collapse All',
                count: '{total}',
                emptyMessage: '<em>No data</em>',
                emptyPanes: 'No SearchPanes',
                loadMessage: 'Loading Search Panes...',
                showMessage: 'Show All',
                title: 'Filters Active - %d'
            },
            layout: 'auto',
            order: [],
            panes: [],
            preSelect: [],
            viewCount: true,
            viewTotal: false
        };
        return SearchPanes;
    }());

    var __extends = (window && window.__extends) || (function () {
        var extendStatics = function (d, b) {
            extendStatics = Object.setPrototypeOf ||
                ({ __proto__: [] } instanceof Array && function (d, b) { d.__proto__ = b; }) ||
                function (d, b) { for (var p in b) if (b.hasOwnProperty(p)) d[p] = b[p]; };
            return extendStatics(d, b);
        };
        return function (d, b) {
            extendStatics(d, b);
            function __() { this.constructor = d; }
            d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
        };
    })();
    var SearchPanesST = /** @class */ (function (_super) {
        __extends(SearchPanesST, _super);
        function SearchPanesST(paneSettings, opts, fromPreInit) {
            if (fromPreInit === void 0) { fromPreInit = false; }
            var _this = this;
            var paneClass;
            if (opts.cascadePanes && opts.viewTotal) {
                paneClass = SearchPaneCascadeViewTotal;
            }
            else if (opts.cascadePanes) {
                paneClass = SearchPaneCascade;
            }
            else if (opts.viewTotal) {
                paneClass = SearchPaneViewTotal;
            }
            _this = _super.call(this, paneSettings, opts, fromPreInit, paneClass) || this;
            var dt = _this.s.dt;
            var loadedFilter = dt.state.loaded();
            var loadFn = function () { return _this._initSelectionListeners(true, loadedFilter && loadedFilter.searchPanes && loadedFilter.searchPanes.selectionList ?
                loadedFilter.searchPanes.selectionList :
                _this.c.preSelect); };
            if (dt.settings()[0]._bInitComplete) {
                loadFn();
            }
            else {
                dt.off('init.dtsps').on('init.dtsps', loadFn);
            }
            return _this;
        }
        /**
         * Ensures that the correct selection listeners are set for selection tracking
         *
         * @param preSelect Any values that are to be preselected
         */
        SearchPanesST.prototype._initSelectionListeners = function (isPreselect, preSelect) {
            if (isPreselect === void 0) { isPreselect = true; }
            if (preSelect === void 0) { preSelect = []; }
            if (isPreselect) {
                this.s.selectionList = preSelect;
            }
            // Set selection listeners for each pane
            for (var _i = 0, _a = this.s.panes; _i < _a.length; _i++) {
                var pane = _a[_i];
                if (pane.s.displayed) {
                    pane.s.dtPane
                        .off('select.dtsp')
                        .on('select.dtsp', this._update(pane))
                        .off('deselect.dtsp')
                        .on('deselect.dtsp', this._updateTimeout(pane));
                }
            }
            // Update on every draw
            this.s.dt.off('draw.dtsps').on('draw.dtsps', this._update());
            // Also update right now as table has just initialised
            this._updateSelectionList();
        };
        /**
         * Retrieve the total values from the server data
         */
        SearchPanesST.prototype._serverTotals = function () {
            for (var _i = 0, _a = this.s.panes; _i < _a.length; _i++) {
                var pane = _a[_i];
                if (pane.s.colOpts.show) {
                    var colTitle = this.s.dt.column(pane.s.index).dataSrc();
                    var blockVT = true;
                    // If any of the counts are not equal to the totals filtering must be active
                    if (this.s.serverData.searchPanes.options[colTitle]) {
                        for (var _b = 0, _c = this.s.serverData.searchPanes.options[colTitle]; _b < _c.length; _b++) {
                            var data = _c[_b];
                            if (data.total !== data.count) {
                                blockVT = false;
                                break;
                            }
                        }
                    }
                    // Set if filtering is present on the pane and populate the data arrays
                    pane.s.filteringActive = !blockVT;
                    pane._serverPopulate(this.s.serverData);
                }
            }
        };
        /**
         * Set's the function that is to be performed when a state is loaded
         *
         * Overrides the method in SearchPanes
         */
        SearchPanesST.prototype._stateLoadListener = function () {
            var _this = this;
            var stateLoadFunction = function (e, settings, data) {
                if (data.searchPanes === undefined) {
                    return;
                }
                // Set the selection list for the panes so that the correct
                // rows can be reselected and in the right order
                _this.s.selectionList =
                    data.searchPanes.selectionList ?
                        data.searchPanes.selectionList :
                        [];
                // Find the panes that match from the state and the actual instance
                if (data.searchPanes.panes) {
                    for (var _i = 0, _a = data.searchPanes.panes; _i < _a.length; _i++) {
                        var loadedPane = _a[_i];
                        for (var _b = 0, _c = _this.s.panes; _b < _c.length; _b++) {
                            var pane = _c[_b];
                            if (loadedPane.id === pane.s.index && pane.s.dtPane) {
                                // Set the value of the searchbox
                                pane.dom.searchBox.val(loadedPane.searchTerm);
                                // Set the value of the order
                                pane.s.dtPane.order(loadedPane.order);
                            }
                        }
                    }
                }
                _this._updateSelectionList();
            };
            this.s.dt.off('stateLoadParams.dtsps', stateLoadFunction).on('stateLoadParams.dtsps', stateLoadFunction);
        };
        /**
         * Remove the function's actions when using cascade
         *
         * Overrides the method in SearchPanes
         */
        SearchPanesST.prototype._updateSelection = function () {
            return;
        };
        /**
         * Returns a function that updates the selection list based on a specific pane
         * Also clears the timeout to stop the deselect from running
         *
         * @param pane the pane that is to have it's selections loaded
         */
        SearchPanesST.prototype._update = function (pane) {
            var _this = this;
            if (pane === void 0) { pane = undefined; }
            return function () {
                if (pane) {
                    clearTimeout(pane.s.deselectTimeout);
                }
                _this._updateSelectionList(pane);
            };
        };
        /**
         * Returns a function that updates the selection list based on a specific pane
         * Also sets a timeout incase a select is about to be made
         *
         * @param pane the pane that is to have it's selections loaded
         */
        SearchPanesST.prototype._updateTimeout = function (pane) {
            var _this = this;
            if (pane === void 0) { pane = undefined; }
            return function () { return pane ?
                pane.s.deselectTimeout = setTimeout(function () { return _this._updateSelectionList(pane); }, 50) :
                _this._updateSelectionList(); };
        };
        /**
         * Updates the selection list to include the latest selections for a given pane
         *
         * @param index The index of the pane that is to be updated
         * @param selected Which rows are selected within the pane
         */
        SearchPanesST.prototype._updateSelectionList = function (paneIn) {
            if (paneIn === void 0) { paneIn = undefined; }
            // Bail if any of these flags are set
            if (this.s.pagingST) {
                // Reset pagingST flag
                this.s.pagingST = false;
                return;
            }
            else if (this.s.updating || paneIn && paneIn.s.serverSelecting) {
                return;
            }
            if (paneIn !== undefined) {
                if (this.s.dt.page.info().serverSide) {
                    paneIn._updateSelection();
                }
                // Get filter values for all of the rows and the selections
                var rows = paneIn.s.dtPane.rows({ selected: true }).data().toArray().map(function (el) { return el.filter; });
                this.s.selectionList = this.s.selectionList.filter(function (selection) { return selection.column !== paneIn.s.index; });
                if (rows.length > 0) {
                    this.s.selectionList.push({
                        column: paneIn.s.index,
                        rows: rows
                    });
                    paneIn.dom.clear.removeClass(this.classes.disabledButton).removeAttr('disabled');
                }
                else {
                    paneIn.dom.clear.addClass(this.classes.disabledButton).attr('disabled', 'true');
                }
                if (this.s.dt.page.info().serverSide) {
                    this.s.dt.draw(false);
                }
            }
            this._remakeSelections();
            this._updateFilterCount();
        };
        /**
         * Remake the selections that were present before new data or calculations have occured
         */
        SearchPanesST.prototype._remakeSelections = function () {
            this.s.updating = true;
            if (!this.s.dt.page.info().serverSide) {
                var tmpSL = this.s.selectionList;
                var anotherFilter = false;
                this.clearSelections();
                this.s.dt.draw(false);
                // When there are no selections present if the length of the data does not match the searched data
                // then another filter is present
                if (this.s.dt.rows().toArray()[0].length > this.s.dt.rows({ search: 'applied' }).toArray()[0].length) {
                    anotherFilter = true;
                }
                this.s.selectionList = tmpSL;
                // Update the rows in each pane
                for (var _i = 0, _a = this.s.panes; _i < _a.length; _i++) {
                    var pane = _a[_i];
                    if (pane.s.displayed) {
                        pane.s.filteringActive = anotherFilter;
                        pane.updateRows();
                    }
                }
                for (var _b = 0, _c = this.s.selectionList; _b < _c.length; _b++) {
                    var selection = _c[_b];
                    var pane = void 0;
                    for (var _d = 0, _e = this.s.panes; _d < _e.length; _d++) {
                        var paneCheck = _e[_d];
                        if (paneCheck.s.index === selection.column) {
                            pane = paneCheck;
                            break;
                        }
                    }
                    if (!pane.s.dtPane) {
                        continue;
                    }
                    var ids = pane.s.dtPane.rows().indexes().toArray();
                    // Select the rows that are present in the selection list
                    for (var i = 0; i < selection.rows.length; i++) {
                        var rowFound = false;
                        for (var _f = 0, ids_1 = ids; _f < ids_1.length; _f++) {
                            var id = ids_1[_f];
                            var currRow = pane.s.dtPane.row(id);
                            var data = currRow.data();
                            if (selection.rows[i] === data.filter) {
                                currRow.select();
                                rowFound = true;
                            }
                        }
                        if (!rowFound) {
                            selection.rows.splice(i, 1);
                            i--;
                        }
                    }
                    pane.s.selections = selection.rows;
                    // If there are no rows selected then don't bother continuing past here
                    // Will just increase processing time and skew the rows that are shown in the table
                    if (selection.rows.length === 0) {
                        continue;
                    }
                    // Update the table to display the current results
                    this.s.dt.draw(false);
                    var filteringActive = false;
                    var filterCount = 0;
                    var prevSelectedPanes = 0;
                    var selectedPanes = 0;
                    // Add the number of all of the filters throughout the panes
                    for (var _g = 0, _h = this.s.panes; _g < _h.length; _g++) {
                        var currPane = _h[_g];
                        if (currPane.s.dtPane) {
                            filterCount += currPane.getPaneCount();
                            if (filterCount > prevSelectedPanes) {
                                selectedPanes++;
                                prevSelectedPanes = filterCount;
                            }
                        }
                    }
                    filteringActive = filterCount > 0;
                    for (var _j = 0, _k = this.s.panes; _j < _k.length; _j++) {
                        var currPane = _k[_j];
                        if (currPane.s.displayed) {
                            // Set the filtering active flag
                            if (anotherFilter || pane.s.index !== currPane.s.index || !filteringActive) {
                                currPane.s.filteringActive = filteringActive || anotherFilter;
                            }
                            else if (selectedPanes === 1) {
                                currPane.s.filteringActive = false;
                            }
                            // Update the rows to show correct counts
                            if (currPane.s.index !== pane.s.index) {
                                currPane.updateRows();
                            }
                        }
                    }
                }
                // Update table to show final search results
                this.s.dt.draw(false);
            }
            else {
                // Identify the last pane to have a change in selection
                var pane = void 0;
                if (this.s.selectionList.length > 0) {
                    pane = this.s.panes[this.s.selectionList[this.s.selectionList.length - 1].column];
                }
                // Update the rows of all of the other panes
                for (var _l = 0, _m = this.s.panes; _l < _m.length; _l++) {
                    var currPane = _m[_l];
                    if (currPane.s.displayed && (!pane || currPane.s.index !== pane.s.index)) {
                        currPane.updateRows();
                    }
                }
            }
            this.s.updating = false;
        };
        return SearchPanesST;
    }(SearchPanes));

    /*! SearchPanes 2.1.2
     * © SpryMedia Ltd - datatables.net/license
     */
    setJQuery$4($);
    setJQuery($);
    setJQuery$3($);
    setJQuery$2($);
    setJQuery$1($);
    var dataTable = $.fn.dataTable;
    // eslint-disable-next-line no-extra-parens
    dataTable.SearchPanes = SearchPanes;
    // eslint-disable-next-line no-extra-parens
    DataTable.SearchPanes = SearchPanes;
    // eslint-disable-next-line no-extra-parens
    dataTable.SearchPanesST = SearchPanesST;
    // eslint-disable-next-line no-extra-parens
    DataTable.SearchPanesST = SearchPanesST;
    // eslint-disable-next-line no-extra-parens
    dataTable.SearchPane = SearchPane;
    // eslint-disable-next-line no-extra-parens
    DataTable.SearchPane = SearchPane;
    // eslint-disable-next-line no-extra-parens
    dataTable.SearchPaneViewTotal = SearchPaneViewTotal;
    // eslint-disable-next-line no-extra-parens
    DataTable.SearchPaneViewTotal = SearchPaneViewTotal;
    // eslint-disable-next-line no-extra-parens
    dataTable.SearchPaneCascade = SearchPaneCascade;
    // eslint-disable-next-line no-extra-parens
    DataTable.SearchPaneCascade = SearchPaneCascade;
    // eslint-disable-next-line no-extra-parens
    dataTable.SearchPaneCascadeViewTotal = SearchPaneCascadeViewTotal;
    // eslint-disable-next-line no-extra-parens
    DataTable.SearchPaneCascadeViewTotal = SearchPaneCascadeViewTotal;
    // eslint-disable-next-line no-extra-parens
    var apiRegister = $.fn.dataTable.Api.register;
    apiRegister('searchPanes()', function () {
        return this;
    });
    apiRegister('searchPanes.clearSelections()', function () {
        return this.iterator('table', function (ctx) {
            if (ctx._searchPanes) {
                ctx._searchPanes.clearSelections();
            }
        });
    });
    apiRegister('searchPanes.rebuildPane()', function (targetIdx, maintainSelections) {
        return this.iterator('table', function (ctx) {
            if (ctx._searchPanes) {
                ctx._searchPanes.rebuild(targetIdx, maintainSelections);
            }
        });
    });
    apiRegister('searchPanes.resizePanes()', function () {
        var ctx = this.context[0];
        return ctx._searchPanes ?
            ctx._searchPanes.resizePanes() :
            null;
    });
    apiRegister('searchPanes.container()', function () {
        var ctx = this.context[0];
        return ctx._searchPanes
            ? ctx._searchPanes.getNode()
            : null;
    });
    DataTable.ext.buttons.searchPanesClear = {
        action: function (e, dt) {
            dt.searchPanes.clearSelections();
        },
        text: 'Clear Panes'
    };
    DataTable.ext.buttons.searchPanes = {
        action: function (e, dt, node, config) {
            var _this = this;
            if (!config._panes) {
                // No SearchPanes on this button yet - initialise and show
                this.processing(true);
                setTimeout(function () {
                    _buttonSourced(dt, node, config);
                    _this.popover(config._panes.getNode(), {
                        align: 'container',
                        span: 'container'
                    });
                    config._panes.rebuild(undefined, true);
                    _this.processing(false);
                }, 10);
            }
            else {
                // Already got SP - show it
                this.popover(config._panes.getNode(), {
                    align: 'container',
                    span: 'container'
                });
                config._panes.rebuild(undefined, true);
            }
        },
        init: function (dt, node, config) {
            dt.button(node).text(config.text || dt.i18n('searchPanes.collapse', 'SearchPanes', 0));
            // If state save is enabled, we need to initialise SP immediately
            // to allow any saved state to be restored. Otherwise we can delay
            // the init until needed by button press
            if (dt.init().stateSave) {
                _buttonSourced(dt, node, config);
            }
        },
        config: {},
        text: ''
    };
    function _buttonSourced(dt, node, config) {
        var buttonOpts = $.extend({
            filterChanged: function (count) {
                dt.button(node).text(dt.i18n('searchPanes.collapse', dt.context[0].oLanguage.searchPanes !== undefined ?
                    dt.context[0].oLanguage.searchPanes.collapse :
                    dt.context[0]._searchPanes.c.i18n.collapse, count));
            }
        }, config.config);
        var panes = buttonOpts && (buttonOpts.cascadePanes || buttonOpts.viewTotal) ?
            new DataTable.SearchPanesST(dt, buttonOpts) :
            new DataTable.SearchPanes(dt, buttonOpts);
        dt.button(node).text(config.text || dt.i18n('searchPanes.collapse', panes.c.i18n.collapse, 0));
        config._panes = panes;
    }
    function _init(settings, options, fromPre) {
        if (options === void 0) { options = null; }
        if (fromPre === void 0) { fromPre = false; }
        var api = new dataTable.Api(settings);
        var opts = options
            ? options
            : api.init().searchPanes || dataTable.defaults.searchPanes;
        var searchPanes = opts && (opts.cascadePanes || opts.viewTotal) ?
            new SearchPanesST(api, opts, fromPre) :
            new SearchPanes(api, opts, fromPre);
        var node = searchPanes.getNode();
        return node;
    }
    // Attach a listener to the document which listens for DataTables initialisation
    // events so we can automatically initialise
    $(document).on('preInit.dt.dtsp', function (e, settings) {
        if (e.namespace !== 'dt') {
            return;
        }
        if (settings.oInit.searchPanes ||
            DataTable.defaults.searchPanes) {
            if (!settings._searchPanes) {
                _init(settings, null, true);
            }
        }
    });
    // DataTables `dom` feature option
    DataTable.ext.feature.push({
        cFeature: 'P',
        fnInit: _init
    });
    // DataTables 2 layout feature
    if (DataTable.ext.features) {
        DataTable.ext.features.register('searchPanes', _init);
    }

})();


return DataTable;
}));


/*! Bootstrap 5 integration for DataTables' SearchPanes
 * © SpryMedia Ltd - datatables.net/license
 */

(function( factory ){
	if ( typeof define === 'function' && define.amd ) {
		// AMD
		define( ['jquery', 'datatables.net-bs5', 'datatables.net-searchpanes'], function ( $ ) {
			return factory( $, window, document );
		} );
	}
	else if ( typeof exports === 'object' ) {
		// CommonJS
		var jq = require('jquery');
		var cjsRequires = function (root, $) {
			if ( ! $.fn.dataTable ) {
				require('datatables.net-bs5')(root, $);
			}

			if ( ! $.fn.dataTable.SearchPanes ) {
				require('datatables.net-searchpanes')(root, $);
			}
		};

		if (typeof window !== 'undefined') {
			module.exports = function (root, $) {
				if ( ! root ) {
					// CommonJS environments without a window global must pass a
					// root. This will give an error otherwise
					root = window;
				}

				if ( ! $ ) {
					$ = jq( root );
				}

				cjsRequires( root, $ );
				return factory( $, root, root.document );
			};
		}
		else {
			cjsRequires( window, jq );
			module.exports = factory( jq, window, window.document );
		}
	}
	else {
		// Browser
		factory( jQuery, window, document );
	}
}(function( $, window, document, undefined ) {
'use strict';
var DataTable = $.fn.dataTable;


$.extend(true, DataTable.SearchPane.classes, {
    buttonGroup: 'btn-group',
    disabledButton: 'disabled',
    narrow: 'col',
    pane: {
        container: 'table'
    },
    paneButton: 'btn btn-light',
    pill: 'badge rounded-pill bg-secondary',
    search: 'form-control search',
    table: 'table table-sm table-borderless',
    topRow: 'dtsp-topRow'
});
$.extend(true, DataTable.SearchPanes.classes, {
    clearAll: 'dtsp-clearAll btn btn-light',
    collapseAll: 'dtsp-collapseAll btn btn-light',
    container: 'dtsp-searchPanes',
    disabledButton: 'disabled',
    panes: 'dtsp-panes dtsp-panesContainer',
    search: DataTable.SearchPane.classes.search,
    showAll: 'dtsp-showAll btn btn-light',
    title: 'dtsp-title',
    titleRow: 'dtsp-titleRow'
});


return DataTable;
}));


