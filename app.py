"""ASTRA: real SPHEREx image discovery and exploratory epoch comparison."""
import io, base64, math
from urllib.parse import urlparse
import requests
import numpy as np
from flask import Flask, request, jsonify
from astropy.io import fits, votable
from astropy.wcs import WCS
from scipy.ndimage import map_coordinates, label, find_objects
from PIL import Image

app = Flask(__name__)
app.config['MAX_CONTENT_LENGTH'] = 15 * 1024 * 1024
SIA = 'https://irsa.ipac.caltech.edu/SIA'
SESSION = requests.Session()
SESSION.headers['User-Agent'] = 'ASTRA-educational-SPHEREx-client/1.0'

def fail(message, code=400): return jsonify(error=message), code

def number(value, lo, hi, name):
    try: n = float(value)
    except (TypeError, ValueError): raise ValueError(f'{name} must be a number')
    if not math.isfinite(n) or not lo <= n <= hi: raise ValueError(f'{name} must be between {lo} and {hi}')
    return n

def allowed_url(url):
    p = urlparse(url)
    return p.scheme == 'https' and p.hostname == 'irsa.ipac.caltech.edu' and p.path.startswith('/ibe/data/spherex/') and p.path.lower().endswith('.fits') and not p.username and not p.password

def rows_from_votable(content):
    table = votable.parse(io.BytesIO(content)).get_first_table().to_table()
    print('IRSA raw rows:', len(table), 'columns:', table.colnames, flush=True)
    out=[]
    for row in table:
        d={}
        for name in table.colnames:
            v=row[name]
            if np.ma.is_masked(v): v=None
            elif isinstance(v, bytes): v=v.decode('utf-8','replace')
            elif isinstance(v, np.generic): v=v.item()
            d[name]=v
        url=str(d.get('access_url') or '')
        if not allowed_url(url): continue
        out.append({'url':url,'title':str(d.get('obs_title') or d.get('obs_id') or 'SPHEREx image'),
                    'date':str(d.get('t_min') or ''),'band_min_m':d.get('em_min'),
                    'band_max_m':d.get('em_max'),'collection':str(d.get('obs_collection') or ''),
                    'id':str(d.get('obs_id') or '')})
    print('IRSA eligible FITS rows:', len(out), flush=True)
    return out

@app.get('/')
def home(): return jsonify(status='online',message='ASTRA API. Open React UI at http://127.0.0.1:5173')

@app.get('/api/spherex/search')
def search():
    try:
        ra=number(request.args.get('ra'),0,360,'RA'); dec=number(request.args.get('dec'),-90,90,'Dec')
        radius=number(request.args.get('radius','0.03'),0.002,0.12,'Radius')
        collection=request.args.get('collection','spherex_qr3')
        if collection not in ('spherex_qr2','spherex_qr3'): return fail('Invalid collection')
        params={'COLLECTION':collection,'POS':f'CIRCLE {ra} {dec} {radius}','RESPONSEFORMAT':'votable','MAXREC':'100'}
        r=SESSION.get(SIA,params=params,timeout=45);r.raise_for_status()
        print('IRSA status:',r.status_code,'type:',r.headers.get('Content-Type'),'preview:',r.text[:160],flush=True)
        results=rows_from_votable(r.content)
        results.sort(key=lambda x:float(x['date']) if x['date'] not in ('','None') else 0)
        return jsonify(results=results[:100],count=len(results),note='Actual IRSA metadata. Choose two observations with overlapping wavelength ranges and distinct times.')
    except ValueError as e: return fail(str(e))
    except Exception as e: return fail('IRSA search failed: '+str(e),502)

def fetch_cutout(url,ra,dec,size):
    if not allowed_url(url): raise ValueError('Only SPHEREx IRSA FITS access URLs are allowed')
    r=SESSION.get(url,params={'center':f'{ra},{dec}','size':str(size),'gzip':'false'},timeout=100,stream=True)
    r.raise_for_status()
    buf=io.BytesIO(); total=0
    for chunk in r.iter_content(131072):
        total+=len(chunk)
        if total>35_000_000: raise ValueError('Cutout exceeds 35 MB; reduce field size')
        buf.write(chunk)
    buf.seek(0)
    with fits.open(buf,memmap=False) as hdus:
        for hdu in hdus:
            if hdu.name.upper()=='IMAGE' and hdu.data is not None and hdu.data.ndim==2:
                return np.array(hdu.data,dtype=np.float32),WCS(hdu.header).celestial
        raise ValueError('No 2D IMAGE extension found in SPHEREx cutout')

def png(arr):
    valid=np.isfinite(arr)
    if not valid.any(): raise ValueError('No valid image pixels')
    lo,hi=np.nanpercentile(arr[valid],[3,99.5]); hi=max(hi,lo+1e-8)
    scaled=np.uint8(np.clip((np.nan_to_num(arr,nan=lo)-lo)/(hi-lo),0,1)**.6*255)
    b=io.BytesIO();Image.fromarray(scaled).save(b,format='PNG')
    return 'data:image/png;base64,'+base64.b64encode(b.getvalue()).decode()

@app.post('/api/spherex/compare')
def compare():
    data=request.get_json(silent=True) or {}
    try:
        ra=number(data.get('ra'),0,360,'RA');dec=number(data.get('dec'),-90,90,'Dec')
        size=number(data.get('size',0.04),0.005,0.12,'Size')
        u1,u2=data.get('before_url',''),data.get('after_url','')
        if u1==u2: return fail('Select two different observations')
        a,wa=fetch_cutout(u1,ra,dec,size); b,wb=fetch_cutout(u2,ra,dec,size)
        # Resample both on a common celestial WCS grid, avoiding false pixel-offset changes.
        h,w=a.shape;h=min(h,700);w=min(w,700);a=a[:h,:w]
        yy,xx=np.indices((h,w),dtype=float)
        lon,lat=wa.pixel_to_world_values(xx,yy)
        bx,by=wb.world_to_pixel_values(lon,lat)
        aligned=map_coordinates(b,[by,bx],order=1,mode='constant',cval=np.nan)
        good=np.isfinite(a)&np.isfinite(aligned)
        if good.sum()<500: raise ValueError('Insufficient overlapping sky pixels; choose another pair')
        # Background/scale matching is only approximate, NOT photometric calibration.
        am=np.nanmedian(a[good]); bm=np.nanmedian(aligned[good]);
        aa=a-am;bb=aligned-bm
        ascale=np.nanpercentile(np.abs(aa[good]),90);bscale=np.nanpercentile(np.abs(bb[good]),90)
        if ascale<1e-12 or bscale<1e-12: raise ValueError('Image has insufficient contrast')
        delta=np.abs(aa/ascale-bb/bscale)
        values=delta[good];threshold=np.percentile(values,99)
        flagged=good&(delta>=threshold)
        overlay=np.zeros((h,w,3),dtype=np.uint8)
        baseline=np.uint8(np.clip((np.nan_to_num(bb,nan=0)/bscale+1)*110,0,255))
        overlay[:]=np.stack([baseline]*3,axis=-1)
        overlay[flagged]=[255,80,85]
        # Connected bright-difference regions: exploratory pixel groups only.
        # These are NOT classified astrophysical objects or confidence scores.
        regions, count = label(flagged)
        candidates = []
        for idx, region in enumerate(find_objects(regions), 1):
            if region is None: continue
            mask = regions[region] == idx
            pixels = int(mask.sum())
            if pixels < 4: continue
            ys, xs = region
            cy = float(ys.start + np.mean(np.where(mask)[0]))
            cx = float(xs.start + np.mean(np.where(mask)[1]))
            try:
                cra, cdec = wa.pixel_to_world_values(cx, cy)
                cra, cdec = float(cra), float(cdec)
            except Exception:
                cra, cdec = None, None
            candidates.append(dict(id=f'ASTRA-REGION-{idx:03d}', classification='Unverified difference region',
                pixel_count=pixels, pixel_x=round(cx, 1), pixel_y=round(cy, 1),
                ra=cra, dec=cdec, bbox=[xs.start, ys.start, xs.stop-xs.start, ys.stop-ys.start]))
        candidates.sort(key=lambda x: x['pixel_count'], reverse=True)
        candidates = candidates[:30]
        out=io.BytesIO();Image.fromarray(overlay).save(out,format='PNG')
        return jsonify(before=png(a),after=png(aligned),overlay='data:image/png;base64,'+base64.b64encode(out.getvalue()).decode(),
                       overlap_pixels=int(good.sum()),flagged_pixels=int(flagged.sum()),
                       threshold=float(threshold), candidates_list=candidates,
                       new_candidates_count=len(candidates), light_curve=[],
                       note='REAL SPHEREx FITS cutouts, WCS-aligned. Red = top 1% normalized differences, NOT verified moving objects or calibrated flux variability.')
    except ValueError as e:return fail(str(e))
    except Exception as e:return fail('Cutout/analysis failed: '+str(e),502)

if __name__=='__main__': app.run(debug=True)
