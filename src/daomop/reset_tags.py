from vos import Client

c = Client()

#for expnum in c.listdir('vos:cfis/solar_system/dbimages', force=True):
for expnum in open('expnum_list.txt').readlines():
    expnum = expnum.strip()
    n = c.get_node('vos:cfis/solar_system/dbimages/{}'.format(expnum), force=True)
    print n.uri
    for prop in n.props:
        if 'stationary' in prop:
           n.props[prop] = None
    c.add_props(n)
    print("Done")
