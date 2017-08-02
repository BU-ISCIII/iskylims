# -*- coding: utf-8 -*-
# Generated by Django 1.11.3 on 2017-07-27 14:11
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('wetlab', '0005_auto_20170727_1332'),
    ]

    operations = [
        migrations.AddField(
            model_name='bioinfo',
            name='department',
            field=models.CharField(default='not assigned', max_length=45),
        ),
        migrations.AddField(
            model_name='bioinfo',
            name='info2',
            field=models.CharField(default='', max_length=45),
        ),
        migrations.AddField(
            model_name='bioinfo',
            name='researcher',
            field=models.CharField(default='', max_length=45),
        ),
        migrations.AddField(
            model_name='bioinfo',
            name='serviceRegistration',
            field=models.CharField(default='not required', max_length=45),
        ),
    ]
